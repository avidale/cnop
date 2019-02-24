

/***
Before execution, we need to specify:
	DGP
	MDLS
	start_iter
	sim_iter
	quiet

For example:
	DGP	= "NOPC"
	MDLS = "OP"
	start_iter	= 1
	sim_iter	= 10*1000
	quiet	= 1
	MIN_CLASS_PERCENTAGE = 0.02
	MIN_CLASS_COUNT = 8
	repeat_dataset = 1

Before execution, we need to call:
	CNOPishModel_definition.ado
	
	gradients.ado
	inflatedOP_estimation_routines.ado
	simulation_routines.do
***/


mata


generate_params_app(DGP, param_true=.,  beta=., a=., gammap=., gamman=., mup=., mun=., ron=., rop=., gamma=., mu=., ro=.)

x1 = st_data(., "spread")
x2 = st_data(., "pb")
x3 = st_data(., "houst")
x4 = st_data(., "gdp")
x	= x1, x2, x3, x4
zp	= x1, x2
zn	= x1, x4
x = J(repeat_dataset, 1, x)
zp = J(repeat_dataset, 1, zp)
zn = J(repeat_dataset, 1, zn)
n = rows(x)
ncat = 5
infcat = 3

/* todo: calculate these things for all rows) */
get_true_me_p_v4(DGP, param_true, x, true_me = ., true_pr = .)

ready = 0
con1 = 0
con2 = 0

for(it = start_iter; it <= 50000; it++){
	/* iteration number, number of good attempts so far, number of converged attempts */
	it, ready, con1, con2
	
	/* Generate data */
	rseed(42+it)

	e0 = rnormal(n,1,0,1)
	e1 = rnormal(n,1,0,1)
	e2 = rnormal(n,1,0,1)
	
	one_simulation(DGP, y=., n, ncat, infcat, x, zp, zn, e0, e1, e2, beta, a, gammap, gamman, mup, mun, ron, rop, gamma, mu, ro)
	
	y_hist = sum(y:==1), sum(y:==2), sum(y:==3), sum(y:==4), sum(y:==5)
	
	if (min(y_hist/rows(y)) < MIN_CLASS_PERCENTAGE  || min(y_hist) < MIN_CLASS_COUNT) {
		y_hist / rows(y)
		"bad data generated, continue another y"
		continue
	}
	
	/* todo: loop over different MDL's and save all parameters into several objects */
	
	estimate_and_get_params_v3(MDLS[1], p1=., s1=., me1=., mese1 = ., pr1 = ., prse1 = ., ll_obs1 = ., acc1 = ., brier1 = ., rps1 = ., aic1 = ., caic1 = ., bic1 = ., lik1 = ., conv1 = ., etime1 = ., eiter1 = ., y=y, x=x, zp=zp, zn=zn, infcat=infcat, quiet=quiet)
	
	estimate_and_get_params_v3(MDLS[2], p2=., s2=., me2=., mese2 = ., pr2 = ., prse2 = ., ll_obs2=. , acc2=., brier2=., rps2=., aic2=., caic2=., bic2=., lik2=., conv2 = ., etime2 = ., eiter2 = ., y=y, x=x, zp=zp, zn=zn, infcat=infcat, quiet=quiet)
	
	if ((conv1 == 0) || (conv2 == 0)) {
		"One of the models did not converge"
		continue
	}
	
	
	pbias1 = pr1 - true_pr 							/* in the end, make it absolute and average */
	prmse1 = (pr1 - true_pr) :^ 2 					/* in the end, average and take root */
	pcovr1 = abs((true_pr - pr1) :/ prse1) :< 1.96 	/*  in the end, average */
	pavgse1 = prse1 								/* in the end, calculate mean predicted se and compare with true */
	pavgm11 = pr1 									/* in the end, calculate true se */
	pavgm21 = pr1 :^ 2 								/* in the end, calculate true se */
	
	mbias1 = me1 - true_me
	mrmse1 = (me1 - true_me) :^ 2
	mcovr1 = abs((true_me - me1) :/ mese1) :< 1.96 
	mavgse1 = mese1 
	mavgm11 = me1 
	mavgm21 = me1 :^ 2
	
	pbias2 = pr2 - true_pr 
	prmse2 = (pr2 - true_pr) :^ 2 
	pcovr2 = abs((true_pr - pr2) :/ prse2) :< 1.96 
	pavgse2 = prse2 
	pavgm12 = pr2 
	pavgm22 = pr2 :^ 2 
	
	mbias2 = me2 - true_me
	mrmse2 = (me2 - true_me) :^ 2
	mcovr2 = abs((true_me - me2) :/ mese2) :< 1.96 
	mavgse2 = mese2 
	mavgm12 = me2 
	mavgm22 = me2 :^ 2
	
	/* same calculation for model parameters */
	if (MDLS[1] == DGP) {
		true_model_p = p1
		true_model_s = s1
	} else if (MDLS[2] == DGP) {
		true_model_p = p2
		true_model_s = s2
	} else {
		/* here we would just throw an error */
	}
	xbias  = true_model_p - param_true'
	xrmse  = (true_model_p - param_true') :^ 2
	xcovr  = abs((true_model_p - param_true') :/ true_model_s) :< 1.96 
	xavgse = true_model_s 
	xavgm1 = true_model_p 
	xavgm2 = true_model_p :^ 2
		
	
	/* likelihood ratio test - we don't conduct it because models are not nested */
	k_1 = cols(p1)
	k_2 = cols(p2)
	/*chi2_pvalue1 = chi2(k_1 - k_2, 2 *(lik1 - lik2))
	chi2_pvalue2 = chi2(k_2 - k_1, 2 *(lik2 - lik2))*/
	
	/* the Vuong part */
	ll_diff = ll_obs1 - ll_obs2
	
	
	mean_diff = mean(ll_diff)
	std_diff = sqrt(variance(ll_diff))
	n_obs = rows(ll_diff)
	vuong = mean_diff / (std_diff / sqrt(n_obs))
	vuongAIC = (mean_diff - (k_1-k_2) / n_obs) / (std_diff / sqrt(n_obs))
	vuongBIC = (mean_diff - (k_1-k_2) * log(n_obs) / (2 * n_obs)) / (std_diff / sqrt(n_obs))
	pvalue = 1-normal(vuong)
	pvalueAIC = 1-normal(vuongAIC)
	pvalueBIC = 1-normal(vuongBIC)
	
	/* sign test */
	n_wins = sum(ll_diff :> 0)
	n_winsAIC = sum((ll_diff :- (k_1-k_2) / n_obs) :> 0)
	n_winsBIC = sum((ll_diff :- (k_1-k_2) * log(n_obs) / (2 * n_obs)) :> 0)
	pr_wins = n_wins / n_obs
	wins_pvalue = binomial(n_obs, n_obs-n_wins, 0.5)
	wins_pvalueAIC = binomial(n_obs, n_obs-n_winsAIC, 0.5)
	wins_pvalueBIC = binomial(n_obs, n_obs-n_winsBIC, 0.5)
	
	crit1 = acc1, brier1, rps1, aic1, caic1, bic1, lik1
	crit2 = acc2, brier2, rps2, aic2, caic2, bic2, lik2
	crit12 = vuong, vuongAIC, vuongBIC, n_wins, n_winsAIC, n_winsBIC, pvalue, pvalueAIC, pvalueBIC, wins_pvalue, wins_pvalueAIC, wins_pvalueBIC
	
	
	cmp1 = acc1 > acc2, brier1 < brier2, rps1 < rps2, aic1 < aic2, caic1 < caic2, bic1 < bic2, lik1 > lik2, pvalue < 0.05, pvalueAIC < 0.05, pvalueBIC < 0.05, wins_pvalue < 0.05, wins_pvalueAIC < 0.05, wins_pvalueBIC < 0.05
	cmp2 = acc2 > acc1, brier2 < brier1, rps2 < rps1, aic2 < aic1, caic2 < caic1, bic2 < bic1, lik2 > lik1, pvalue > 0.95, pvalueAIC > 0.95, pvalueBIC > 0.95, wins_pvalue > 0.95, wins_pvalueAIC > 0.95, wins_pvalueBIC > 0.95
	
	if (ready == 0) {
		all_cmp1 = cmp1
		all_cmp2 = cmp2
		
		all_crit1  = crit1
		all_crit2  = crit2
		all_crit12 = crit12
		
		all_pbias1 = pbias1 
		all_prmse1 = prmse1
		all_pcovr1 = pcovr1
		all_pavgse1 = pavgse1
		all_pavgm11 = pavgm11
		all_pavgm21 = pavgm21
		
		all_mbias1 = mbias1 
		all_mrmse1 = mrmse1
		all_mcovr1 = mcovr1
		all_mavgse1 = mavgse1
		all_mavgm11 = mavgm11
		all_mavgm21 = mavgm21
		
		all_pbias2 = pbias2 
		all_prmse2 = prmse2
		all_pcovr2 = pcovr2
		all_pavgse2 = pavgse2
		all_pavgm12 = pavgm12
		all_pavgm22 = pavgm22
		
		all_mbias2 = mbias2
		all_mrmse2 = mrmse2
		all_mcovr2 = mcovr2
		all_mavgse2 = mavgse2
		all_mavgm12 = mavgm12
		all_mavgm22 = mavgm22
		
		all_xbias = xbias
		all_xrmse = xrmse
		all_xcovr = xcovr
		all_xavgse = xavgse
		all_xavgm1 = xavgm1
		all_xavgm2 = xavgm2
		
		all_y = y
		
	} else {
		all_cmp1 = all_cmp1 \ cmp1
		all_cmp2 = all_cmp2 \ cmp2
		
		all_crit1  = all_crit1 \ crit1
		all_crit2  = all_crit2 \ crit2
		all_crit12 = all_crit12 \ crit12
		
		all_pbias1 = all_pbias1 \ pbias1 
		all_prmse1 = all_prmse1 \ prmse1
		all_pcovr1 = all_pcovr1 \ pcovr1
		all_pavgse1 = all_pavgse1 \ pavgse1
		all_pavgm11 = all_pavgm11 \ pavgm11
		all_pavgm21 = all_pavgm21 \ pavgm21
		
		all_mbias1 = all_mbias1 \ mbias1 
		all_mrmse1 = all_mrmse1 \ mrmse1
		all_mcovr1 = all_mcovr1 \ mcovr1
		all_mavgse1 = all_mavgse1 \ mavgse1
		all_mavgm11 = all_mavgm11 \ mavgm11
		all_mavgm21 = all_mavgm21 \ mavgm21
		
		all_pbias2 = all_pbias2 \ pbias2
		all_prmse2 = all_prmse2 \ prmse2
		all_pcovr2 = all_pcovr2 \ pcovr2
		all_pavgse2 = all_pavgse2 \ pavgse2
		all_pavgm12 = all_pavgm12 \ pavgm12
		all_pavgm22 = all_pavgm22 \ pavgm22
		
		all_mbias2 = all_mbias2 \ mbias2 
		all_mrmse2 = all_mrmse2 \ mrmse2
		all_mcovr2 = all_mcovr2 \ mcovr2
		all_mavgse2 = all_mavgse2 \ mavgse2
		all_mavgm12 = all_mavgm12 \ mavgm12
		all_mavgm22 = all_mavgm22 \ mavgm22
		
		all_xbias = all_xbias \ xbias
		all_xrmse = all_xrmse \ xrmse
		all_xcovr = all_xcovr \ xcovr
		all_xavgse = all_xavgse \ xavgse
		all_xavgm1 = all_xavgm1 \ xavgm1
		all_xavgm2 = all_xavgm2 \ xavgm2
		
		all_y = all_y \ y
		
	}
	
	
	ready = ready + 1
	con1 = con1 + conv1
	con2 = con2 + conv2
	if(ready >= sim_iter){
		break
	}
	
	/*
	if (conv == 1) {
		"another successful row added!"
		if (ready == 0) {
			allpa = p
			allse = s
			allme = me
			allms = mese
			allpr = pr
			allps = prse
			allco = conv
			allet = etime
			allyh = y_hist
			allei = eiter
			allit = it
		} else {
			allpa = allpa \ p
			allse = allse \ s
			allme = allme \ me
			allms = allms \ mese
			allpr = allpr \ pr
			allps = allps \ prse
			allco = allco \ conv
			allet = allet \ etime
			allyh = allyh \ y_hist
			allei = allei \ eiter
			allit = allit \ it
		}
	} else {
		"not converged, continue another y"
	}
	*/
}

it
rows(sim_iter)


true_categories = all_y /* just repeat y multiple times */

effp1_new = calucalate_metrics(all_pbias1, all_prmse1, all_pcovr1, all_pavgse1, all_pavgm11, all_pavgm21, 1, true_categories, 1)
effm1_new = calucalate_metrics(all_mbias1, all_mrmse1, all_mcovr1, all_mavgse1, all_mavgm11, all_mavgm21, 1, true_categories, 4)

effp2_new = calucalate_metrics(all_pbias2, all_prmse2, all_pcovr2, all_pavgse2, all_pavgm12, all_pavgm22, 1, true_categories, 1)
effm2_new = calucalate_metrics(all_mbias2, all_mrmse2, all_mcovr2, all_mavgse2, all_mavgm12, all_mavgm22, 1, true_categories, 4)

/* for marginal effects, reshape every row into 4 rows, by variable */
effm1_new = effm1_new[,(1,2,3,4,5,6)] \ effm1_new[,(7,8,9,10,11,12)] \ effm1_new[,(13,14,15,16,17,18)] \ effm1_new[,(19,20,21,22,23,24)]
effm2_new = effm2_new[,(1,2,3,4,5,6)] \ effm2_new[,(7,8,9,10,11,12)] \ effm2_new[,(13,14,15,16,17,18)] \ effm2_new[,(19,20,21,22,23,24)]
/* now put each variable into its own part */
effm1_new = effm1_new[(1,5,9,13),] \ effm1_new[(2,6,10,14),] \ effm1_new[(3,7,11,15),] \ effm1_new[(4,8,12,16),]
effm2_new = effm2_new[(1,5,9,13),] \ effm2_new[(2,6,10,14),] \ effm2_new[(3,7,11,15),] \ effm2_new[(4,8,12,16),]

effx_new = calucalate_metrics(all_xbias, all_xrmse, all_xcovr, all_xavgse, all_xavgm1, all_xavgm2, 0)


prmse1 = (colsum(addRepeatedSelectColumns(all_pavgm11 - getDummies(all_y), true_categories, 1) :^2) / rows(all_y)) :^ 0.5
prmse2 = (colsum(addRepeatedSelectColumns(all_pavgm21 - getDummies(all_y), true_categories, 1) :^2) / rows(all_y)) :^ 0.5


/* this is the basic table

column1 = effp1, effm1, colsum(all_cmp1) / rows(all_cmp1)

column2 = effp2, effm2, colsum(all_cmp2) / rows(all_cmp2)

column1', column2'
 */

/* this is the mean values of one-model scores (first part of comparison) */

colsum(all_crit1)' / rows(all_crit1), colsum(all_crit2)' / rows(all_crit2)

/* this is the mean values of two-model scores and their p-values (second part of comparison) */

colsum(all_crit12)' / rows(all_crit12)

/*
sim = allpa, allse, allme, allms, allpr, allps, allco, allet, allyh, allei, allit

sim = allpa, allse, allme, allms, allpr, allps, allco, allet, allyh
*/

/* 
fname = "mc_" + strlower(DGP) + "_with_" + strlower(MDL) + "_" + strofreal(n) + "_partial_"+strofreal(start_iter)+"to"+strofreal(sim_iter)+".matamatrix"
fh = fopen(fname, "w")
fputmatrix(fh, sim)
fclose(fh)

*/


end



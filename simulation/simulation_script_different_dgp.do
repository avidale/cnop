

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

Before execution, we need to call:
	CNOPishModel_definition.ado
	gradients.ado
	inflatedOP_estimation_routines.ado
	simulation_routines.do
***/


mata


generate_params_app(DGP, param_true=.,  beta=., a=., gammap=., gamman=., mup=., mun=., ron=., rop=., gamma=., mu=., ro=.)

ready = 0
con1 = 0
con2 = 0

for(it = start_iter; it <= 50000; it++){
	/* iteration number, number of good attempts, number of converged attempts */
	it, ready, con1, con2
	
	/* Generate data */
	it2 = 0
	rseed(42+it)
	
	x1 = st_data(., "spread")
	x2 = st_data(., "pb")
	x3 = st_data(., "houst")
	x4 = st_data(., "gdp")
	n = rows(x4)

	e0 = rnormal(n,1,0,1)
	e1 = rnormal(n,1,0,1)
	e2 = rnormal(n,1,0,1)
	
	x	= x1, x2, x3, x4
	zp	= x1, x2
	zn	= x1, x4

	ncat = 5
	infcat = 3
	
	one_simulation(DGP, y=., n, ncat, infcat, x, zp, zn, e0, e1, e2, beta, a, gammap, gamman, mup, mun, ron, rop, gamma, mu, ro)
	
	y_hist = sum(y:==1), sum(y:==2), sum(y:==3), sum(y:==4), sum(y:==5)
	
	if (min(y_hist/rows(y))<0.03) {
		y_hist / rows(y)
		"bad data generated, continue another y"
		continue
	}
	
	/* todo: loop over different MDL's and save all parameters into several objects */
	
	estimate_and_get_params_v3(MDLS[1], p1=., s1=., me1=., mese1 = ., pr1 = ., prse1 = ., ll_obs1 = ., acc1 = ., brier1 = ., rps1 = ., aic1 = ., caic1 = ., bic1 = ., lik1 = ., conv1 = ., etime1 = ., eiter1 = ., y=y, x=x, zp=zp, zn=zn, infcat=infcat, quiet=quiet)
	
	estimate_and_get_params_v3(MDLS[2], p2=., s2=., me2=., mese2 = ., pr2 = ., prse2 = ., ll_obs2=. , acc2=., brier2=., rps2=., aic2=., caic2=., bic2=., lik2=., conv2 = ., etime2 = ., eiter2 = ., y=y, x=x, zp=zp, zn=zn, infcat=infcat, quiet=quiet)
	
	
	get_true_me_p_v3(DGP, param_true, true_me = ., true_pr = .)
	
	pbias1 = pr1 - true_pr /* todo: in the end, make it absolute and average */
	prmse1 = (pr1 - true_pr) :^ 2 /* todo: in the end, average and take root */
	pcovr1 = abs((true_pr - pr1) :/ prse1) :< 1.96 /* todo: in the end, average */
	pavgse1 = prse1 /* todo: in the end, calculate mean predicted se and compare with true */
	pavgm11 = pr1 /* todo: in the end, calculate true se */
	pavgm21 = pr1 :^ 2 /* todo: in the end, calculate true se */
	
	/* todo: 
	
	for each model:
		calculate bias, RMSE, CR and bias of s.e. for choice probabilities
		calculate bias, RMSE, CR and bias of s.e. for their marginal effects
	
	compare and mark the winner:
		# of correct predictions
		# brier score
		# ranked prob score
		# AIC
		# CAIC
		# BIC
		# Likelihood
		
		# Vuong
		# Vuong AIC
		# Vuong BIC
		
		# Sign
		# Sign AIC
		# Sign BIC
	*/
	
	/* the Vuong part */
	ll_diff = ll_obs1 - ll_obs2
	k_1 = cols(p1)
	k_2 = cols(p2)
	
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
	
	
	n1 = acc1 > acc2, brier1 < brier2, rps1 < rps2, aic1 < aic2, caic1 < caic2, bic1 < bic2, lik1 > lik2, vuong > 0, vuongAIC > 0, vuongBIC > 0, n_wins > n_obs/2, n_winsAIC > n_obs/2, n_winsBIC > n_obs/2
	n2 = acc2 > acc1, brier2 < brier1, rps2 < rps1, aic2 < aic1, caic2 < caic1, bic2 < bic1, lik2 > lik1, vuong < 0, vuongAIC < 0, vuongBIC < 0, n_wins < n_obs/2, n_winsAIC < n_obs/2, n_winsBIC < n_obs/2
	n1
	
	if (ready == 0) {
		nn1 = n1
		nn2 = n2
	} else {
		nn1 = nn1 \ n1
		nn2 = nn2 \ n2
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

/* this row may be incorrect because params may have different shape
param_true,  mean(allpa)'
*/
/*
param_true
mean(allpa)'

mean(allco)
*/

colsum(nn1) / rows(nn1)

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



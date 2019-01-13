

mata


function generate_params_app(dgp, param_true, | beta, a, gammap, gamman, mup, mun, ron, rop, gamma, mu, ro) {
	/* todo: choose best parameters for both processes */
	/*
	Variables are spread pb houst gdp, target is rate_change
	*/
	if (dgp == "OP") {
		beta = (1.5742 \ 0.9262 \ 1.3732 \ 0.2391)
		a = (0.4659 \ 1.8382 \  4.8360 \ 6.3318)
		param_true = beta \ a
	} else if (dgp == "CNOP") {
		/* this was the full model, and we want a restricted one 
		beta 	= ( 2.0730 \ 1.3727 \ 5.7837 \ 0.4287)
		a =       (10.4269 \ 13.2016)
		gammap 	= ( 3.0316 \ 3.6308 \-3.9376 \ 0.6806)
		gamman 	= ( 1.0466 \ 0.5979 \-0.5597 \ 0.1381)
		mup =     (-4.7426 \ 3.4933)
		mun =     (-1.9382 \-0.5704)
		*/
		
		beta 	= ( 2.1063 \ 1.6285 \ 5.3114 \ 0.3810) /* spread, pb, houst, gdp */
		a =       ( 9.1035 \12.3481)
		gammap 	= ( 1.8097 \ 2.6201) /* spread, pb */
		gamman 	= ( 1.0729 \ 0.1777) /* spread, gdp */
		mup =     (-1.4818 \ 3.5091) 
		mun =     (-0.6374 \ 0.7570)

		param_true = beta \ a \ gammap \ mup \ gamman \ mun 
	} else {
		"DGP " + dgp + " is not supported for application-like parameter generation"
		param_true = .
	}
}


function generate_params(dgp, param_true, | beta, a, gammap, gamman, mup, mun, ron, rop, gamma, mu, ro) {
	if (dgp == "NOP"){
		beta 	= (0.6 \ 0.4) 	
		a = (0.21 \ 2.19)
		gammap 	= (0.3 \ 0.9)
		gamman 	= (0.2 \ 0.3)
		mup = (0.68)
		mun = (-0.17)

		param_true = beta \ a \ gammap \ mup \ gamman \ mun 
	}
	if (dgp == "NOPC"){
		beta 	= (0.6 \ 0.4) 	
		a = (0.21 \ 2.19)
		gammap 	= (0.3 \ 0.9)
		gamman 	= (0.2 \ 0.3)
		mup = ( 1.31)
		mun = (-0.5)
		rop = 0.6
		ron = 0.3

		param_true = beta \ a \ gammap \ mup \ gamman \ mun \ ron \ rop
	}
	if (dgp == "MIOPR"){
		beta 	= (0.6 \ 0.8) 	
		a = (0.45)
		gamma 	= (0.5 \ 0.6)
		mu = (-1.45 \ -0.55 \ 0.75 \ 1.65)

		param_true = beta \ a \ gamma \ mu 
	}
	if (dgp == "MIOPRC"){
		beta 	= (0.6 \ 0.8) 	
		a = (0.45)
		gamma 	= (0.5 \ 0.6)
		mu = (-1.18 \ -0.33 \ 0.9 \ 1.76)
		ro = 0.5

		param_true = beta \ a \ gamma \ mu \ ro
	}
	if (dgp == "CNOP"){
		
		beta 	= (0.6 \ 0.4) 	
		a = (0.9 \ 1.5)
		gammap 	= (0.3 \ 0.9)
		gamman 	= (0.2 \ 0.3)
		mup = ( 0.02 \ 1.28)
		mun = (-0.67 \ 0.36)

		param_true = beta \ a \ gammap \ mup \ gamman \ mun 
	}
	if (dgp == "CNOPC"){
		beta 	= (0.6 \ 0.4) 	
		a = (0.9 \ 1.5)
		gammap 	= (0.3 \ 0.9)
		gamman 	= (0.2 \ 0.3)
		mup = ( 0.49 \ 1.67)
		mun = (-0.88 \ 0.12)
		rop = 0.6
		ron = 0.3

		param_true = beta \ a \ gammap \ mup \ gamman \ mun \ ron \ rop
	}
}

function one_simulation(dgp, y, n, ncat, infcat, x, zp, zn, e0, e1, e2, beta, a, gammap, gamman, mup, mun, ron, rop, gamma, mu, ro) {
	if (dgp == "OP"){
		int0 = x  * beta + e0
		y = J(n, 1, 1)
		for( i = 1; i <= ncat-1; i++){
			y = y + (int0 :> a[i])
		}
	} 
	else if (dgp == "NOP"){
		eps0 = e0
		epsp = e1
		epsn = e2
		
		int0 = x  * beta + eps0
		intp = zp * gammap + epsp
		intn = zn * gamman + epsn

		rp = (int0 :> a[2])
		rn = (int0 :< a[1])
		r0 = 1:-rp:-rn

		y = J(n, 1, 1) + (rp+r0):*(infcat-1) + rp * 1
		for( i = 1; i <= infcat-2; i++){
			y = y + (intn :> mun[i]):*rn
		}
		for( i = 1; i <= ncat-infcat-1; i++){
			y = y + (intp :> mup[i]):*rp
		}
	} 
	else if (dgp == "NOPC"){
		eps0 = e0
		epsp = e1 * sqrt(1 - rop^2) + eps0 * rop
		epsn = e2 * sqrt(1 - ron^2) + eps0 * ron

		int0 = x  * beta + eps0
		intp = zp * gammap + epsp
		intn = zn * gamman + epsn

		rp = (int0 :> a[2])
		rn = (int0 :< a[1])
		r0 = 1:-rp:-rn

		y = J(n, 1, 1) + (rp+r0):*(infcat-1) + rp * 1
		for( i = 1; i <= infcat-2; i++){
			y = y + (intn :> mun[i]):*rn
		}
		for( i = 1; i <= ncat-infcat-1; i++){
			y = y + (intp :> mup[i]):*rp
		}
	}
	else if (dgp == "MIOPR"){
		eps0 = e0
		eps1 = e1

		int0 = x  * beta + eps0
		int1 = zp * gamma + eps1

		r1 = (int0 :> a[1])
		r0 = 1:-r1

		y = J(n, 1, 1) + r0 * (infcat-1)
		for( i = 1; i <= ncat-1; i++){
			y = y + (int1 :> mu[i]):*r1
		}
	}
	else if (dgp == "MIOPRC"){
		eps0 = e0
		eps1 = e1 * sqrt(1 - ro^2) + eps0 * ro

		int0 = x  * beta + eps0
		int1 = zp * gamma + eps1

		r1 = (int0 :> a[1])
		r0 = 1:-r1

		y = J(n, 1, 1) + r0 * (infcat-1)
		for( i = 1; i <= ncat-1; i++){
			y = y + (int1 :> mu[i]):*r1
		}
	}
	else if (dgp == "CNOP"){
		eps0 = e0
		epsp = e1
		epsn = e2

		int0 = x  * beta + eps0
		intp = zp * gammap + epsp
		intn = zn * gamman + epsn

		rp = (int0 :> a[2])
		rn = (int0 :< a[1])
		r0 = 1:-rp:-rn

		y	= J(n, 1, 1)
		y = y + (rp+r0):*(infcat-1)
		for( i = 1; i <= infcat-1; i++){
			y = y + (intn :> mun[i]):*rn
		}
		for( i = 1; i <= ncat-infcat; i++){
			y = y + (intp :> mup[i]):*rp
		}
	}
	else if (dgp == "CNOPC"){
		eps0 = e0
		epsp = e1 * sqrt(1 - rop^2) + eps0 * rop
		epsn = e2 * sqrt(1 - ron^2) + eps0 * ron

		int0 = x  * beta + eps0
		intp = zp * gammap + epsp
		intn = zn * gamman + epsn

		rp = (int0 :> a[2])
		rn = (int0 :< a[1])
		r0 = 1:-rp:-rn

		y	= J(n, 1, 1)
		y = y + (rp+r0):*(infcat-1)
		for( i = 1; i <= infcat-1; i++){
			y = y + (intn :> mun[i]):*rn
		}
		for( i = 1; i <= ncat-infcat; i++){
			y = y + (intp :> mup[i]):*rp
		}
	} 
	else {
		"DGP " + dgp + " is not supported for y generation"
	}
}


function estimate_and_get_params_v2(dgp, p, s, me, mese, pr, prse, conv, etime, eiter, y, x, zp, zn, infcat, quiet) {
	class CNOPModel scalar mod
	if (dgp == "OP") {
		mod = estimateOP(y, x, quiet) 
	} 
	else if (dgp == "NOP") {
		mod = estimateNOP(y, x, zp, zn, infcat, quiet) 
	} 
	else if (dgp == "NOPC") {
		mod = estimateNOPC(y, x, zp, zn, infcat, quiet) 
	} 
	else if (dgp == "MIOPR") {
		mod = estimateMIOPR(y, x, zp, infcat, quiet) 
	} 
	else if (dgp == "MIOPRC") {
		mod = estimateMIOPRC(y, x, zp, infcat, quiet) 
	} 
	else if (dgp == "CNOP") {
		mod = estimateCNOP(y, x, zp, zn, infcat, quiet) 
	} 
	else if (dgp == "CNOPC") {
		mod = estimateCNOPC(y, x, zp, zn, infcat, quiet) 
	} else {
		"Don't know how to estimate model specified as " + dgp
	}
	p = mod.params'
	s = mod.se'
	mod.corresp = (1,1,0 \ 0,1,1 \ 1,0,1)
	conv = mod.converged
	etime = mod.etime
	eiter = mod.iterations
	if (conv == 1) {
		/* todo: decide whether we need the argument  (0,0,1) which has been the third */
		me_se = generalMEwithSE((2,0,0), mod, 1)
		pr_se = generalPredictWithSE((2,0,0),mod, 1)
	} else {
		me_se = J(6,5,.)
		pr_se = J(2,5,.)
	}
	me = rowshape(me_se[(1,2,3),], 1)
	mese = rowshape(me_se[(4,5,6),], 1)
	pr = pr_se[1,]
	prse = pr_se[2,]
	
	
}


function estimate_and_get_params_v3(dgp, p, s, me, mese, pr, prse, ll_obs, acc, briers, rps, aic, caic, bic, lik, conv, etime, eiter, y, x, zp, zn, infcat, quiet) {
	class CNOPModel scalar mod
	if (dgp == "OP") {
		mod = estimateOP(y, x, quiet) 
	} 
	else if (dgp == "NOP") {
		mod = estimateNOP(y, x, zp, zn, infcat, quiet) 
	} 
	else if (dgp == "NOPC") {
		mod = estimateNOPC(y, x, zp, zn, infcat, quiet) 
	} 
	else if (dgp == "MIOPR") {
		mod = estimateMIOPR(y, x, zp, infcat, quiet) 
	} 
	else if (dgp == "MIOPRC") {
		mod = estimateMIOPRC(y, x, zp, infcat, quiet) 
	} 
	else if (dgp == "CNOP") {
		mod = estimateCNOP(y, x, zp, zn, infcat, quiet) 
	} 
	else if (dgp == "CNOPC") {
		mod = estimateCNOPC(y, x, zp, zn, infcat, quiet) 
	} else {
		"Don't know how to estimate model specified as " + dgp
	}
	p = mod.params'
	s = mod.se'
	mod.corresp = (1,1,1,1 \ 1,1,0,0 \ 1,0,0,1)
	conv = mod.converged
	etime = mod.etime
	eiter = mod.iterations
	
	/*
	column median is (0.063, 0, 1.47, 4.9)
	
	*/
	
	pvt_x = (0.063, 0, 1.47, 4.9)
	
	if (conv == 1) {
		/* todo: decide whether we need the argument  (0,0,1) which has been the third */
		me_se = generalMEwithSE(pvt_x, mod, 1)
		pr_se = generalPredictWithSE(pvt_x, mod, 1)
	} else {
		me_se = J(6,5,.)
		pr_se = J(2,5,.)
	}
	me = rowshape(me_se[(1,2,3,4),], 1)
	mese = rowshape(me_se[(5,6,7,8),], 1)
	pr = pr_se[1,]
	prse = pr_se[2,]
	
	ll_obs = mod.ll_obs
	acc    = mod.accuracy
	briers = mod.brier_score
	rps	   = mod.ranked_probability_score
	aic    = mod.AIC
	bic    = mod.BIC
	caic   = mod.CAIC
	lik    = mod.logLik
}




function get_null_params(dgp, p, s, me, mese, pr, prse, conv, etime, y, x, zp, zn, infcat, quiet) {	
	if (dgp == "NOP"){
		p = J(1, 10, .)
	}
	if (dgp == "NOPC"){
		p = J(1, 12, .)
	}
	if (dgp == "MIOPR"){
		p = J(1, 10, .)
	}
	if (dgp == "MIOPRC"){
		p = J(1, 11, .)
	}
	if (dgp == "CNOP"){
		p = J(1, 12, .)
	}
	if (dgp == "CNOPC"){
		p = J(1, 14, .)
	}
	s = p
	conv = 0
	etime = 0
	
	me_se = J(6,5,.)
	pr_se = J(2,5,.)
	me = rowshape(me_se[(1,2,3),], 1)
	mese = rowshape(me_se[(4,5,6),], 1)
	pr = pr_se[1,]
	prse = pr_se[2,]
	
}

function get_null_params_v2(dgp, p, s, me, mese, pr, prse, conv, etime, eiter, y, x, zp, zn, infcat, quiet) {	
	if (dgp == "NOP"){
		p = J(1, 10, .)
	}
	if (dgp == "NOPC"){
		p = J(1, 12, .)
	}
	if (dgp == "MIOPR"){
		p = J(1, 10, .)
	}
	if (dgp == "MIOPRC"){
		p = J(1, 11, .)
	}
	if (dgp == "CNOP"){
		p = J(1, 12, .)
	}
	if (dgp == "CNOPC"){
		p = J(1, 14, .)
	}
	s = p
	conv = 0
	etime = 0
	eiter = 0
	
	me_se = J(6,5,.)
	pr_se = J(2,5,.)
	me = rowshape(me_se[(1,2,3),], 1)
	mese = rowshape(me_se[(4,5,6),], 1)
	pr = pr_se[1,]
	prse = pr_se[2,]
	
}

function get_true_me_p(dgp, par_true, _returnedME, _returnedP){
	
	class CNOPModel scalar true_model
	true_model.model_class 	= dgp
	true_model.params 	= par_true
	
	true_model.ncat 	= 5
	true_model.infcat 	= 3
	true_model.corresp	= (1,1,0 \ 0,1,1 \ 1,0,1)
	
	/* todo: decide whether we needed the 4th argument (0,0,1) */
	generalME(par_true', (2,0,0), true_model, 1, _returnedME = .)
	
	_returnedP = generalPredict((2,0,0), true_model, 1)
}

function get_true_me_p_v3(dgp, par_true, _returnedME, _returnedP){
	
	class CNOPModel scalar true_model
	true_model.model_class 	= dgp
	true_model.params 	= par_true
	
	true_model.ncat 	= 5
	true_model.infcat 	= 3
	true_model.corresp	= (1,1,1,1 \ 1,1,0,0 \ 1,0,0,1)
	pvt_x = (0.063, 0, 1.47, 4.9)
	
	generalME(par_true', pvt_x, true_model, 1, _returnedME = .)
	
	_returnedP = generalPredict(pvt_x, true_model, 1)
}


function colMedians(matrix x) {
	k = cols(x)
	n = rows(x)
	results = J(1, k, 0)
	for (i = 1; i <= k; i++) {
		y = x[,i]
		_sort(y,1)
		if (mod(n,2) == 1) {
			results[i] = y[(n+1)/2]
		} else {
			results[i] = 0.5*(y[n/2]+y[n/2+1])
		}
	}
	return(results)
}

	
/*
	if (dgp == "NOP"){
	}
	if (dgp == "NOPC"){
	}
	if (dgp == "MIOPR"){
	}
	if (dgp == "MIOPRC"){
	}
	if (dgp == "CNOP"){
	}
	if (dgp == "CNOPC"){
	}
*/

end







mata
/*

generate_params("NOPC", param_true=.,  beta=., a=., gammap=., gamman=., mup=., mun=., ron=., rop=., gamma=., mu=., ro=.)

rseed(42)
n = 1000


x1	= rnormal(n,1,0,1) :+ 2
x2	= rnormal(n,1,0,1)
x3	= runiform(n,1)
x3	= (x3 :> 0.70) :- (x3 :< 0.30)

e0 = rnormal(n,1,0,1)
e1 = rnormal(n,1,0,1)
e2 = rnormal(n,1,0,1)

x	= x1, x2
zp	= x2, x3
zn	= x1, x3

ncat = 5
infcat = 3

one_simulation("NOPC", y=., n, ncat, infcat, x, zp, zn, e0, e1, e2, beta, a, gammap, gamman, mup, mun, ron, rop, gamma, mu, ro)

estimate_and_get_params("NOPC", p=., s=., me=., mese = ., pr = ., prse = ., conv = ., etime = ., y=y, x=x, zp=zp, zn=zn, infcat=infcat, quiet=0)

*/
end

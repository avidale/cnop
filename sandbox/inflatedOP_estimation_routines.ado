version 12

// requires CNOPishModel_definition.ado
// requires gradients.ado

mata

function matrix_mse(residual) {
	result = mean(rowsum(residual :* residual))
	return(result)
}

function running_rowsum(input_matrix) {
	result = input_matrix
	for(i=2; i<=cols(input_matrix); i++) {
		result[,i] = result[,i-1] + result[,i]
	}
	return(result)
}

// calculate and remember various postestimation statistics
class CNOPModel scalar describeModel(class CNOPModel scalar model, params, covMat, covMat_rob, maxLik, n, q, prob_obs) {
	
	model.params = params
	model.se		= sqrt(diagonal(covMat))
	model.t			= abs(params :/ model.se)
	model.se_rob	= sqrt(diagonal(covMat_rob))
	model.t_rob		= abs(params :/ model.se_rob)
	
	model.AIC	= -2 * maxLik + 2 * rows(params) 
	model.BIC	= -2 * maxLik + ln(n) * rows(params)
	model.CAIC	= -2 * maxLik + (1 + ln(n)) * rows(params)
	model.AICc	= model.AIC + 2 * rows(params) * (rows(params) + 1) / (n - rows(params) - 1)
	model.HQIC	= -2 * maxLik + 2*rows(params)*ln(ln(n))
	model.logLik0 	= sum(log(q :* mean(q)))
	model.R2 	= 1 - maxLik /  model.logLik0
	
	model.brier_score = matrix_mse(prob_obs - q)
	model.ranked_probability_score = matrix_mse(running_rowsum(prob_obs) - running_rowsum(q))
	
	model.V	= covMat
	model.V_rob	= covMat_rob
	model.logLik	= maxLik
	model.probabilities = prob_obs
	model.ll_obs = log(rowsum(prob_obs :* q))
	
	return (model)
}


// workfunction that returns OP model (as an instance of CNOPModel class)
class CNOPModel scalar estimateOP(y, x, |quiet, startvalues, xbar, dummies, robust, who){
	// incomplete: change argument list
	/*
	feed q instead of y ?
	feed_vsop
	*/
	
	starttime = clock(c("current_time"),"hms")
	
	if (rows(quiet) < 1) {
		quiet = 0
	}
	// identify model dimensions
	
	n	= rows(x)
	k	= cols(x)
	allcat	= uniqrows(y)
	ncat	= rows(allcat)
	
	// complete missing values
	if( (cols(xbar) == 0) | (rows(xbar) < k)){
		xbar	= colsum(x) :/ n
	}
	if( (cols(dummies) == 0) | (length(dummies) < k)){
		dummies	= J(1, k, 0)
	}
	if(cols(robust) == 0){
		robust	= 0
	}
	
	
	
	// compute categories
	q = J(n, ncat, 0)
	for(i=1; i<=ncat; i++){
		q[.,i]	=(y :== allcat[i])
	}
	
	// initial values
	start_param	= startvalues
	if(rows(start_param) <= 1){
		start_mu	= invnormal(runningsum(mean(q))[1::ncat-1])';
		start_b		= invsym(x'*x)*x'*y;
		start_param	= start_b \ start_mu 
	}
	
	// optimize
	// incomplete: set optimization parameters
	S = optimize_init()
	if(quiet){
		optimize_init_tracelevel(S , "none")
		optimize_init_verbose(S, 0)
	}
	optimize_init_argument(S, 1, x)
	optimize_init_argument(S, 2, q)
	optimize_init_argument(S, 3, ncat)
	optimize_init_evaluator(S, &_op_optim())
	optimize_init_evaluatortype(S, "gf1") // added!!
	optimize_init_params(S, (start_param'))
	
	if (cols(who) > 0 && who != .) {
		optimize_init_cluster(S, who) // added!!!
	}
	errorcode = _optimize(S)
	if (errorcode == 0) {
		params = optimize_result_params(S)'
	} else {
		"dd: OP final optimization encountered error code " + strofreal(errorcode)
		params = optimize_result_params(S)'
	}
	
	// extract optimization results
	maxLik	= optimize_result_value(S)
	grad 	= optimize_result_gradient(S)
	covMat	= optimize_result_V_oim(S)	// added!
	retCode	= optimize_result_errortext(S)
	
	
	// calculate robust variance
	covMat_rob	= optimize_result_V_robust(S)		// in contrast with the source, here I get it automatically ??? no, the variable "who" matters !!! 
	 
	g	= Jacop(params, x, q, ncat) // gradient for all observations: n \times rows(params)
	/*
	ss	= J(rows(params), rows(params), 0)
	for (i=1; i<=max(who);i++){
		sel = select(g, who :== i)
		ss  = ss + sel' * sel; 
	}
	covMat_rob =covMat * ss * covMat; 
	*/
	
	// calculate model statistics
	se		= sqrt(diagonal(covMat))
	tstat	= abs(params :/ se)
	se_rob		= sqrt(diagonal(covMat_rob))
	tstat_rob	= abs(params :/ se_rob)

	// calculate predicted probabilities
	pred_prob	= MLop(params, x, q, ncat, 1)
	
	// calculate predicted option ?????

	// calculate ME
	me		= op_me(params, xbar, q, ncat)
	
	// calculate SE for ME
	mese	= J(k * ncat,1, 0) // to be reshaped
	mese_rob	= J(k * ncat,1, 0) // to be reshaped
	
	// I want to find derivative of all ME's with respect to all parameters
	D = deriv_init()
	deriv_init_evaluator(D, &_op_me_deriv())
	deriv_init_evaluatortype(D, "t")
	deriv_init_argument(D, 1, xbar)
	deriv_init_argument(D, 2, q)
	deriv_init_argument(D, 3, ncat)
	deriv_init_params(D, params') // a row vector
	dydx = deriv(D, 1)
	// in columns - params
	// in rows - elements of rowshape(me[,1::cols(x)],1) - i.e. k * ncat effects

	for(i = 1; i<=rows(dydx); i++){
		mese[i]	= dydx[i,] * covMat * dydx[i,]'
		mese_rob[i]	= dydx[i,] * covMat_rob * dydx[i,]'
	}
	mese	= colshape(mese, ncat)
	mese_rob	= colshape(mese_rob, ncat)
	
	// calculate t-statistics for ME
	met		= abs(me) :/ sqrt(mese)
	met_rob	= abs(me) :/ sqrt(mese_rob)
	
	// just for convenience
	gama = params[1::k]
	mu = params[(k+1)::(k+ncat-1)]
	
	
	// pack results
	class CNOPModel scalar model 
	model.model_class = "OP"
	
	model.n	= n
	model.k	= k
	model.ncat	= ncat
	model.gamma	= gama
	model.mu	= mu
	model.allcat = allcat
	
	model.params	= params
	model.se		= se
	model.t			= tstat
	model.se_rob	= se_rob
	model.t_rob		= tstat_rob
	
	model.me		= me
	model.met		= met
	model.mese		= mese
	
	model.AIC	= -2 * maxLik + 2 * rows(params) 
	model.BIC	= -2 * maxLik + ln(n) * rows(params)
	model.CAIC	= -2 * maxLik + (1 + ln(n)) * rows(params)
	model.AICc	= model.AIC + 2 * rows(params) * (rows(params) + 1) / (n - rows(params) - 1)
	model.HQIC	= -2 * maxLik + 2*rows(params)*ln(ln(n))
	model.logLik0 	= sum(log(q :* mean(q)))
	model.R2 	= 1 - maxLik /  model.logLik0
	
	model.V	= covMat
	model.V_rob	= covMat_rob
	model.logLik	= maxLik
	model.probabilities	= pred_prob
	
	model.retCode	= retCode
	model.etime		= -1 // to be determined via Stata functions ???
	
	model.retCode = optimize_result_errortext(S)
	model.error_code = errorcode
	model.etime = clock(c("current_time"),"hms") - starttime
	model.converged = optimize_result_converged(S)
	
	return(model)
}



// workfunction that returns NOP model (as an instance of CNOPModel class)
class CNOPModel scalar estimateNOP(y, x, zp, zn, infcat, |quiet, startvalues, robust, who, corresp, lambda, maxiter, ptol, vtol, nrtol){
	if (rows(quiet) < 1) {
		quiet = 0
	}
	if (args() < 11 || rows(lambda) < 1 || lambda ==.) {
		external CNOP_GLOBAL_CONST_LAMBDA
		lambda = CNOP_GLOBAL_CONST_LAMBDA
	}
	if (args() < 12 || rows(maxiter) < 1 || maxiter ==.) {
		external CNOP_GLOBAL_CONST_MAXITER
		maxiter = CNOP_GLOBAL_CONST_MAXITER
	}
	if (args() < 13 || rows(ptol) < 1 || ptol ==.) {
		external CNOP_GLOBAL_CONST_PTOL
		ptol = CNOP_GLOBAL_CONST_PTOL
	}
	if (args() < 14 || rows(vtol) < 1 || vtol ==.) {
		external CNOP_GLOBAL_CONST_VTOL
		vtol = CNOP_GLOBAL_CONST_VTOL
	}
	if (args() < 15 || rows(nrtol) < 1 || nrtol ==.) {
		external CNOP_GLOBAL_CONST_NRTOL
		nrtol = CNOP_GLOBAL_CONST_NRTOL
	}
	
	starttime = clock(c("current_time"),"hms")
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	allcat = uniqrows(y)
	ncat = rows(allcat)
	ncatp = sum(allcat :> infcat)
	ncatn = sum(allcat :< infcat)
	infcat_index = selectindex(allcat :== infcat)
	
	q = J(n, ncat, 0)
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== allcat[i])
	}
	q0 = (y:==infcat)
	qp = (y:>infcat)
	qn = (y:<infcat)
	q3 = q0 , qp, qn
	
	start_param	= startvalues //  avoid overwriting external argument
	
	parlen =  (kx+2+kzp+ncatp+kzn+ncatn - 2)
	if (rows(start_param) < parlen && rows(start_param) > 0 && start_param != .) {
		"Vector of initial values must have length "+ strofreal(parlen)
		start_param = .
	}
	if(rows(start_param) < 1 || start_param == .) {
		if(!quiet){
			"Upper-level decision"
		}
		params1 = coeffOP(x, q3, 3, quiet, ., lambda, maxiter, ptol, vtol, nrtol)
		
		fltr = (y:>infcat)
		if(!quiet){
			"Lower-level decision for y>0"
		}
		y2 = select(y, fltr)
		q2 = select(q, fltr)
		q2 = select(q2, (J(1,ncatn+1,0),J(1,ncatp,1))) 
		zp2 = select(zp, fltr)
		
		params2 = coeffOP(zp2, q2, ncatp, quiet, ., lambda, maxiter, ptol, vtol, nrtol)
		if(!quiet){
			"Lower-level decision for y<0"
		}
		fltr = (y:<infcat)
		y2 = select(y, fltr)
		q2 = select(q, fltr)
		q2 = select(q2, (J(1,ncatn,1),J(1,ncatp+1,0))) 
		zn2 = select(zn, fltr)
		
		params3 = coeffOP(zn2, q2, ncatn, quiet, ., lambda, maxiter, ptol, vtol, nrtol)
		
		start_param = params1 \ params2 \ params3
	}
	if(!quiet){
		"Estimating NOP with exogenous switching"
	}
	
	// initially transform arguments to avoid inequality constraints
	_nop_params(start_param, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = .)
	coded_param = b \ codeIncreasingSequence(a) \ gp \ codeIncreasingSequence(mup) \ gn \ codeIncreasingSequence(mun)
	/* if coded params contain mistake, replace them by total zero */
	if (max(coded_param :==.)  > 0){
		coded_param = J(rows(coded_param), cols(coded_param), 0)
		if(!quiet){
			"Initial parameters contain error. Replaced with zero"
		}
	}
	
	// make a few attempts to estimate parameters;
	estimation_successful = 0
	for (estimation_attempt = 0; estimation_attempt <= 7; estimation_attempt++) {
		if(!quiet) {
			"Attempt " + strofreal(estimation_attempt)
		}
		if (estimation_attempt == 0) {
			initial_coded_param = coded_param
			opt_method = "nr"
		}
		if (estimation_attempt == 1) {
			initial_coded_param = coded_param
			opt_method = "bhhh"
		}
		if (estimation_attempt == 2) {
			initial_coded_param = coded_param
			opt_method = "dfp"
		}
		if (estimation_attempt == 3) {
			initial_coded_param = coded_param
			opt_method = "bfgs"
		}
		if (estimation_attempt == 4) {
			initial_coded_param = coded_param * 0
			opt_method = "nr"
		}
		if (estimation_attempt == 5) {
			initial_coded_param = coded_param * 0
			opt_method = "bhhh"
		}
		if (estimation_attempt == 6) {
			initial_coded_param = coded_param * 0
			opt_method = "dfp"
		}
		if (estimation_attempt == 7) {
			initial_coded_param = coded_param * 0
			opt_method = "bfgs"
		}
		singularHmethod = "hybrid" 
		
	
		S = optimize_init()
		if(quiet){
			optimize_init_tracelevel(S , "none")
			optimize_init_verbose(S, 0)
		}
		optimize_init_argument(S, 1, x)
		optimize_init_argument(S, 2, zp)
		optimize_init_argument(S, 3, zn)
		optimize_init_argument(S, 4, q)
		optimize_init_argument(S, 5, ncat)
		optimize_init_argument(S, 6, infcat_index)
		optimize_init_argument(S, 7, 1) // coded
		optimize_init_evaluator(S, &_nop_optim())
		optimize_init_evaluatortype(S, "gf0") // not concave if use GF1
		optimize_init_conv_maxiter(S, maxiter)
		optimize_init_params(S, initial_coded_param')
		optimize_init_conv_ptol(S, ptol)
		optimize_init_conv_vtol(S, vtol)
		optimize_init_conv_nrtol(S, nrtol)
		optimize_init_singularHmethod(S, singularHmethod)
		optimize_init_conv_warning(S, "off") 
		optimize_init_technique(S, opt_method)
		if (cols(who) > 0 && who != .) {
			optimize_init_cluster(S, who) // added!!!
		}
		errorcode = _optimize(S)
		retCode	= optimize_result_errortext(S)
		convg	= optimize_result_converged(S)
		params = optimize_result_params(S)'
		iterations = optimize_result_iterations(S)
		
		if (errorcode == 0) {
			if (convg == 1) {
				estimation_successful = 1
				break; // if optimization is successful, then exit
			} else {
				if(!quiet) {
					"Convergence not achieved"
				}
			}
		} else {
			if(!quiet) {
				"NOP final optimization encountered error code " + strofreal(errorcode) + ": " + retCode
			}
		}
	}
	if (estimation_successful == 0) {
		"Sorry, but despite all attempts, estimation of CNOP parameters did not converge. Turn verbosity on to see the details."
		"Perhaps, there are too few data for such a complex model."
		"Error code is " + strofreal(errorcode) + ": " + retCode
		"Convergence status is " + strofreal(convg)
	}
		
	// After optimization, transform argument back
	_nop_params(params, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = .)
	params = b \ decodeIncreasingSequence(a) \ gp \ decodeIncreasingSequence(mup) \ gn \ decodeIncreasingSequence(mun)
	
	S2 = optimize_init()
	optimize_init_evaluator(S2, &_nop_optim())
	optimize_init_argument(S2, 1, x)
	optimize_init_argument(S2, 2, zp)
	optimize_init_argument(S2, 3, zn)
	optimize_init_argument(S2, 4, q)
	optimize_init_argument(S2, 5, ncat)
	optimize_init_argument(S2, 6, infcat_index)
	optimize_init_argument(S2, 7, 0)
	optimize_init_params(S2, params')
	
	optimize_init_evaluatortype(S2, "gf0")
	optimize_init_conv_maxiter(S2, maxiter)
	optimize_init_conv_ptol(S2, ptol)
	optimize_init_conv_vtol(S2, vtol)
	optimize_init_conv_nrtol(S2, nrtol)
	optimize_init_singularHmethod(S2, "hybrid")
	optimize_init_conv_warning(S2, "off")
	optimize_init_technique(S2, "nr") 
	if (cols(who) > 0 && who != .) {
		optimize_init_cluster(S2, who) // added!!!
	}
	
	errorcode2 = _optimize_evaluate(S2)
	
	
	if(!quiet){
		"Estimation completed"
		""
	}
	
	if (estimation_successful == 0) {
		/* When estimation is not successful, robust covatiance matrix cannot be calculated, and ordinary covariance matrix is . */
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = covMat
	} else {
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = optimize_result_V_robust(S2)
	}
	// calculate predicted probabilities for each observation
	prob_obs = MLnop(params, x, zp, zn, q, ncat, infcat_index, 1)
	
	class CNOPModel scalar model 
	model.model_class = "NOP"
	
	model.n	= n
	model.k	= kx
	model.ncat	= ncat
	model.ncatn = ncatn
	model.ncatp = ncatp
	model.infcat = infcat_index
	model.allcat = allcat
	
	model = describeModel(model, params, covMat, covMat_rob, maxLik, n, q, prob_obs)
	
	model.retCode = retCode
	model.error_code = errorcode
	model.etime = clock(c("current_time"),"hms") - starttime
	model.converged = convg
	model.iterations = iterations
	
	return(model)
	
}







// workfunction that returns ZIOP model (as an instance of CNOPModel class)
class CNOPModel scalar estimateMIOPR(y, x, z, infcat, |quiet, startvalues, robust, who, corresp, lambda, maxiter, ptol, vtol, nrtol) {
	// categories of y must be coded as 1, 2, ... ncat, so that infcat is index of value and of position at the same time
	
	if (rows(quiet) < 1) {
		quiet = 0
	}
	if (args() < 10 || rows(lambda) < 1 || lambda ==.) {
		external CNOP_GLOBAL_CONST_LAMBDA
		lambda = CNOP_GLOBAL_CONST_LAMBDA
	}
	if (args() < 11 || rows(maxiter) < 1 || maxiter ==.) {
		external CNOP_GLOBAL_CONST_MAXITER
		maxiter = CNOP_GLOBAL_CONST_MAXITER
	}
	if (args() < 12 || rows(ptol) < 1 || ptol ==.) {
		external CNOP_GLOBAL_CONST_PTOL
		ptol = CNOP_GLOBAL_CONST_PTOL
	}
	if (args() < 13 || rows(vtol) < 1 || vtol ==.) {
		external CNOP_GLOBAL_CONST_VTOL
		vtol = CNOP_GLOBAL_CONST_VTOL
	}
	if (args() < 14 || rows(nrtol) < 1 || nrtol ==.) {
		external CNOP_GLOBAL_CONST_NRTOL
		nrtol = CNOP_GLOBAL_CONST_NRTOL
	}
	
	starttime = clock(c("current_time"),"hms")
	n	= rows(x)
	kx	= cols(x)
	kz	= cols(z)
	allcat = uniqrows(y)
	ncat = rows(allcat)
	infcat_index = selectindex(allcat :== infcat)
	
	// compute categories
	q = J(n, ncat, 0)
	for(i=1; i<=ncat; i++) {
		q[.,i] = (y :== allcat[i])
	}
	
	// initial values
	start_param	= startvalues //  avoid overwriting external argument
	
	parlen = (kx+kz + ncat)
	if (rows(start_param) < parlen && rows(start_param) > 0 && start_param != .) {
		"Vector of initial values must have length "+ strofreal(parlen)
		start_param = .
	}
	if(rows(start_param) < 1 || start_param == .) {
		if(!quiet){
			"Computing starting values:"
			"Independent regime equation"
		}
		q0	=	(y :== infcat)
		params1 = coeffOP(x, (q0, 1 :- q0), 2, quiet, ., lambda, maxiter, ptol, vtol, nrtol)
		if(!quiet){
			"Independent outcome equation"
		}
		params2 = coeffOP(z, q, ncat, quiet, ., lambda, maxiter, ptol, vtol, nrtol)
			
		start_param	= params1 \ params2
	}
	
	if(!quiet){
		"Estimating ZIOP-2 with exogenous switching"
	}
	
	// initially transform arguments to avoid inequality constraints
	_miopr_params(start_param, kx, kz, ncat, b = ., a = ., g = ., mu = .)
	coded_param = b \ codeIncreasingSequence(a) \ g \ codeIncreasingSequence(mu) 
	/* if coded params contain mistake, replace them by total zero */
	if (max(coded_param :==.)  > 0){
		coded_param = J(rows(coded_param), cols(coded_param), 0)
		if(!quiet){
			"Initial parameters contain error. Replaced with zero"
		}
	}
	
	// make a few attempts to estimate parameters;
	estimation_successful = 0
	for (estimation_attempt = 0; estimation_attempt <= 7; estimation_attempt++) {
		if(!quiet) {
			"Attempt " + strofreal(estimation_attempt)
		}
		if (estimation_attempt == 0) {
			initial_coded_param = coded_param
			opt_method = "nr"
		}
		if (estimation_attempt == 1) {
			initial_coded_param = coded_param
			opt_method = "bhhh"
		}
		if (estimation_attempt == 2) {
			initial_coded_param = coded_param
			opt_method = "dfp"
		}
		if (estimation_attempt == 3) {
			initial_coded_param = coded_param
			opt_method = "bfgs"
		}
		if (estimation_attempt == 4) {
			initial_coded_param = coded_param * 0
			opt_method = "nr"
		}
		if (estimation_attempt == 5) {
			initial_coded_param = coded_param * 0
			opt_method = "bhhh"
		}
		if (estimation_attempt == 6) {
			initial_coded_param = coded_param * 0
			opt_method = "dfp"
		}
		if (estimation_attempt == 7) {
			initial_coded_param = coded_param * 0
			opt_method = "bfgs"
		}
		singularHmethod = "hybrid" 
	
	
		S = optimize_init()
		if(quiet){
			optimize_init_tracelevel(S , "none")
			optimize_init_verbose(S, 0)
		}
		optimize_init_argument(S, 1, x)
		optimize_init_argument(S, 2, z)
		optimize_init_argument(S, 3, q)
		optimize_init_argument(S, 4, ncat)
		optimize_init_argument(S, 5, infcat_index)
		optimize_init_argument(S, 6, 1) // coded
		optimize_init_evaluator(S, &_miopr_optim())
		
		optimize_init_evaluatortype(S, "gf0") // not concave if use GF1
		optimize_init_conv_maxiter(S, maxiter)
		optimize_init_params(S, initial_coded_param')
		optimize_init_conv_ptol(S, ptol)
		optimize_init_conv_vtol(S, vtol)
		optimize_init_conv_nrtol(S, nrtol)
		optimize_init_singularHmethod(S, singularHmethod)
		optimize_init_conv_warning(S, "off") 
		optimize_init_technique(S, opt_method)
		if (cols(who) > 0 && who != .) {
			optimize_init_cluster(S, who) // added!!!
		}
		errorcode = _optimize(S)
		retCode	= optimize_result_errortext(S)
		convg	= optimize_result_converged(S)
		params = optimize_result_params(S)'
		iterations = optimize_result_iterations(S)
		
		if (errorcode == 0) {
			if (convg == 1) {
				estimation_successful = 1
				break; // if optimization is successful, then exit
			} else {
				if(!quiet) {
					"Convergence not achieved"
				}
			}
		} else {
			if(!quiet) {
				"ZIOP-2 final optimization encountered error code " + strofreal(errorcode) + ": " + retCode
			}
		}
	}
	
	if (estimation_successful == 0) {
		"Sorry, but despite all attempts, estimation of ZIOP parameters did not converge. Turn verbosity on to see the details."
		"Perhaps, there are too few data for such a complex model."
		"Error code is " + strofreal(errorcode) + ": " + retCode
		"Convergence status is " + strofreal(convg)
	}
	
	// After optimization, transform argument back
	_miopr_params(params, kx, kz, ncat, b = ., a = ., g = ., mu = .)
	params = b \ decodeIncreasingSequence(a) \ g \ decodeIncreasingSequence(mu) 
	
	S2 = optimize_init()
	optimize_init_argument(S2, 1, x)
	optimize_init_argument(S2, 2, z)
	optimize_init_argument(S2, 3, q)
	optimize_init_argument(S2, 4, ncat)
	optimize_init_argument(S2, 5, infcat_index)
	optimize_init_argument(S2, 6, 0) // flag that params are coded to avoid inequality constraints
	optimize_init_evaluator(S2, &_miopr_optim())
	optimize_init_evaluatortype(S2, "gf0")
	optimize_init_conv_maxiter(S2, maxiter)
	if (cols(who) > 0 && who != .) {
		optimize_init_cluster(S2, who) 
	}
	
	optimize_init_conv_ptol(S2, ptol)
	optimize_init_conv_vtol(S2, vtol)
	optimize_init_conv_nrtol(S2, nrtol)
	optimize_init_singularHmethod(S2, "hybrid")
	optimize_init_conv_warning(S2, "off") // show that convergence not achieved
	optimize_init_technique(S2, "nr") 
	
	optimize_init_params(S2, params')
	
	errorcode2 = _optimize_evaluate(S2)
	
	if(!quiet){
		"Estimation completed"
		""
	}
	//"TEST: CALCULATION REALLY COMPLETED"
	//estimation_successful
	//params
	if (estimation_successful == 0) {
		/* When estimation is not successful, robust covatiance matrix cannot be calculated, and ordinary covariance matrix is . */
		
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = covMat
	} else {
		//"TEST: MAXLIK IS"
		maxLik	= optimize_result_value(S2)
		//maxLik
		//"TEST: GRAD IS"
		grad 	= optimize_result_gradient(S2)
		//grad
		//"TEST: COVMAT IS"
		covMat	= optimize_result_V(S2)
		//covMat
		//"TEST: ROB IS"
		covMat_rob = optimize_result_V_robust(S2)
		//covMat_rob
	}
	// calculate predicted probabilities for each observation
	prob_obs = MLmiopr(params, x, z, q, ncat, infcat_index, 1)
	
	class CNOPModel scalar model 
	model.model_class = "MIOPR"
	
	model.n	= n
	model.k	= kx + kz
	model.ncat	= ncat
	model.infcat = infcat_index
	model.allcat = allcat
	
	model = describeModel(model, params, covMat, covMat_rob, maxLik, n, q, prob_obs)
	
	model.retCode = retCode
	model.error_code = errorcode
	model.etime = clock(c("current_time"),"hms") - starttime
	model.converged = convg
	model.iterations = iterations
	
	//"TEST: EXIT MIOPR ESTIMATION"
	
	return(model)
}

// workfunction that returns ZIOP-correlated model (as an instance of CNOPModel class)
class CNOPModel scalar estimateMIOPRC(y, x, z, infcat, |quiet, startvalues, robust, who, corresp, lambda, maxiter, ptol, vtol, nrtol) {
	// categories of y must be coded as 1, 2, ... ncat, so that infcat is index of value and of position at the same time
	
	if (rows(quiet) < 1) {
		quiet = 0
	}
	if (args() < 10 || rows(lambda) < 1 || lambda ==.) {
		external CNOP_GLOBAL_CONST_LAMBDA
		lambda = CNOP_GLOBAL_CONST_LAMBDA
	}
	if (args() < 11 || rows(maxiter) < 1 || maxiter ==.) {
		external CNOP_GLOBAL_CONST_MAXITER
		maxiter = CNOP_GLOBAL_CONST_MAXITER
	}
	if (args() < 12 || rows(ptol) < 1 || ptol ==.) {
		external CNOP_GLOBAL_CONST_PTOL
		ptol = CNOP_GLOBAL_CONST_PTOL
	}
	if (args() < 13 || rows(vtol) < 1 || vtol ==.) {
		external CNOP_GLOBAL_CONST_VTOL
		vtol = CNOP_GLOBAL_CONST_VTOL
	}
	if (args() < 14 || rows(nrtol) < 1 || nrtol ==.) {
		external CNOP_GLOBAL_CONST_NRTOL
		nrtol = CNOP_GLOBAL_CONST_NRTOL
	}
	
	starttime = clock(c("current_time"),"hms")
	n	= rows(x)
	kx	= cols(x)
	kz	= cols(z)
	allcat = uniqrows(y)
	ncat = rows(allcat)
	infcat_index = selectindex(allcat :== infcat)
	
	q = J(n, ncat, 0)
	for(i=1; i<=ncat; i++) {
		q[.,i] = (y :== allcat[i])
	}
	
	start_param	= startvalues //  avoid overwriting external argument
	
	parlen = (kx+kz+ncat+1)
	if (rows(start_param) < parlen && rows(start_param) > 0 && start_param != .) {
		"Vector of initial values must have length "+ strofreal(parlen)
		start_param = .
	}
	if(rows(start_param) < 1 || start_param == .) {
		if(!quiet) {
			"Run estimation of two-part zero-inflated ordered probit (ZIOP-2) model with endogenous switching"
		}
		class CNOPModel scalar initial_model 
		initial_model = estimateMIOPR(y, x, z, infcat, quiet, ., ., ., ., lambda, maxiter, ptol, vtol, nrtol)
		if(!quiet) {
			"Estimation of ZIOP-2 with exogenous switching successful"
		}
		start_param = initial_model.params \ 0
		/* Here I could have started looping through different ro, but I just set ro = 0, and it's okay; more quick and stable. */
	}
	
	if(!quiet){
		"Estimating ZIOP-2 with endogenous switching"
	}
	
	// initially transform arguments to avoid inequality constraints
	_mioprc_params(start_param, kx, kz, ncat, b = ., a = ., g = ., mu = ., ro = .)
	
	coded_param = b \ codeIncreasingSequence(a) \ g \ codeIncreasingSequence(mu) \ logit((ro+1)/2)
	/* if coded params contain mistake, replace them by total zero */
	if (max(coded_param :==.)  > 0){
		coded_param = J(rows(coded_param), cols(coded_param), 0)
		if(!quiet){
			"Initial parameters contain error. Replaced with zero"
		}
	}
	
	// make a few attempts to estimate parameters;
	estimation_successful = 0
	for (estimation_attempt = 0; estimation_attempt <= 7; estimation_attempt++) {
		if(!quiet) {
			"Attempt " + strofreal(estimation_attempt)
		}
		if (estimation_attempt == 0) {
			initial_coded_param = coded_param
			opt_method = "nr"
		}
		if (estimation_attempt == 1) {
			initial_coded_param = coded_param
			opt_method = "bhhh"
		}
		if (estimation_attempt == 2) {
			initial_coded_param = coded_param
			opt_method = "dfp"
		}
		if (estimation_attempt == 3) {
			initial_coded_param = coded_param
			opt_method = "bfgs"
		}
		if (estimation_attempt == 4) {
			initial_coded_param = coded_param * 0
			opt_method = "nr"
		}
		if (estimation_attempt == 5) {
			initial_coded_param = coded_param * 0
			opt_method = "bhhh"
		}
		if (estimation_attempt == 6) {
			initial_coded_param = coded_param * 0
			opt_method = "dfp"
		}
		if (estimation_attempt == 7) {
			initial_coded_param = coded_param * 0
			opt_method = "bfgs"
		}
		singularHmethod = "hybrid" 
		
		S = optimize_init()
		if(quiet){
			optimize_init_tracelevel(S , "none")
			optimize_init_verbose(S, 0)
		}
		optimize_init_argument(S, 1, x)
		optimize_init_argument(S, 2, z)
		optimize_init_argument(S, 3, q)
		optimize_init_argument(S, 4, ncat)
		optimize_init_argument(S, 5, infcat_index)
		optimize_init_argument(S, 6, 1) // flag that params are coded to avoid inequality constraints
		optimize_init_evaluator(S, &_mioprc_optim())
		optimize_init_evaluatortype(S, "gf0") // not concave if use GF1
		optimize_init_conv_maxiter(S, maxiter)
		optimize_init_params(S, initial_coded_param')
		optimize_init_conv_ptol(S, ptol)
		optimize_init_conv_vtol(S, vtol)
		optimize_init_conv_nrtol(S, nrtol)
		optimize_init_singularHmethod(S, singularHmethod)
		optimize_init_conv_warning(S, "off") 
		optimize_init_technique(S, opt_method)
		if (cols(who) > 0 && who != .) {
			optimize_init_cluster(S, who) // added!!!
		}
		errorcode = _optimize(S)
		retCode	= optimize_result_errortext(S)
		convg	= optimize_result_converged(S)
		params = optimize_result_params(S)'
		iterations = optimize_result_iterations(S)
		
		if (errorcode == 0) {
			if (convg == 1) {
				estimation_successful = 1
				break; // if optimization is successful, then exit
			} else {
				if(!quiet) {
					"Convergence not achieved"
				}
			}
		} else {
			if(!quiet) {
				"ZIOP-2 final optimization encountered error code " + strofreal(errorcode) + ": " + retCode
			}
		}
	}
	if (estimation_successful == 0) {
		"Sorry, but despite all attempts, estimation of ZIOPC parameters did not converge. Turn verbosity on to see the details."
		"Perhaps, there are too few data for such a complex model."
		"Error code is " + strofreal(errorcode) + ": " + retCode
		"Convergence status is " + strofreal(convg)
	}
	
	// After optimization, transform argument back
	_mioprc_params(params, kx, kz, ncat, b = ., a = ., g = ., mu = ., ro = .)
	params = b \ decodeIncreasingSequence(a) \ g \ decodeIncreasingSequence(mu) \ invlogit(ro) * 2-1
	
	S2 = optimize_init()
	optimize_init_argument(S2, 1, x)
	optimize_init_argument(S2, 2, z)
	optimize_init_argument(S2, 3, q)
	optimize_init_argument(S2, 4, ncat)
	optimize_init_argument(S2, 5, infcat_index)
	optimize_init_argument(S2, 6, 0) // flag that params are coded to avoid inequality constraints
	optimize_init_evaluator(S2, &_mioprc_optim())
	optimize_init_evaluatortype(S2, "gf0")
	optimize_init_conv_maxiter(S2, maxiter)
	if (cols(who) > 0 && who != .) {
		optimize_init_cluster(S2, who) // added!!!
	}
	optimize_init_conv_ptol(S2, ptol)
	optimize_init_conv_vtol(S2, vtol)
	optimize_init_conv_nrtol(S2, nrtol)
	optimize_init_singularHmethod(S2, "hybrid")
	optimize_init_conv_warning(S2, "off") // show that convergence not achieved
	optimize_init_technique(S2, "nr") 
	
	optimize_init_params(S2, params')
	
	errorcode = _optimize_evaluate(S2)
	
	if(!quiet){
		"Estimation completed"
		""
	}
	
	if (estimation_successful == 0) {
		/* When estimation is not successful, robust covatiance matrix cannot be calculated, and ordinary covariance matrix is . */
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = covMat
	} else {
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = optimize_result_V_robust(S2)
	}
	// calculate predicted probabilities for each observation
	prob_obs = MLmioprc(params, x, z, q, ncat, infcat_index, 1)
	
	class CNOPModel scalar model 
	model.model_class = "MIOPRC"
	
	model.n	= n
	model.k	= kx + kz
	model.ncat	= ncat
	model.infcat = infcat_index
	model.allcat = allcat
	
	model = describeModel(model, params, covMat, covMat_rob, maxLik, n, q, prob_obs)
	
	model.retCode = retCode
	model.error_code = errorcode
	model.etime = clock(c("current_time"),"hms") - starttime
	model.converged = convg
	model.iterations = iterations
	
	return(model)
}


// workfunction that returns NOP-correlated model (as an instance of CNOPModel class)
class CNOPModel scalar estimateNOPC(y, x, zp, zn, infcat, |quiet, startvalues, robust, who, corresp, lambda, maxiter, ptol, vtol, nrtol){
	if (rows(quiet) < 1) {
		quiet = 0
	}
	if (args() < 11 || rows(lambda) < 1 || lambda ==.) {
		external CNOP_GLOBAL_CONST_LAMBDA
		lambda = CNOP_GLOBAL_CONST_LAMBDA
	}
	if (args() < 12 || rows(maxiter) < 1 || maxiter ==.) {
		external CNOP_GLOBAL_CONST_MAXITER
		maxiter = CNOP_GLOBAL_CONST_MAXITER
	}
	if (args() < 13 || rows(ptol) < 1 || ptol ==.) {
		external CNOP_GLOBAL_CONST_PTOL
		ptol = CNOP_GLOBAL_CONST_PTOL
	}
	if (args() < 14 || rows(vtol) < 1 || vtol ==.) {
		external CNOP_GLOBAL_CONST_VTOL
		vtol = CNOP_GLOBAL_CONST_VTOL
	}
	if (args() < 15 || rows(nrtol) < 1 || nrtol ==.) {
		external CNOP_GLOBAL_CONST_NRTOL
		nrtol = CNOP_GLOBAL_CONST_NRTOL
	}
	
	
	starttime = clock(c("current_time"),"hms")
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	allcat = uniqrows(y)
	ncat = rows(allcat)
	ncatp = sum(allcat :> infcat)
	ncatn = sum(allcat :< infcat)
	infcat_index = selectindex(allcat :== infcat)
	
	q = J(n, ncat, 0)
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== allcat[i])
	}
	q0 = (y:==infcat)
	qp = (y:>infcat)
	qn = (y:<infcat)
	q3 = q0 , qp, qn
	
	start_param	= startvalues //  avoid overwriting external argument
	
	parlen =  (kx+2+kzp+ncatp+kzn+ncatn)
	if (rows(start_param) < parlen && rows(start_param) > 0 && start_param != .) {
		"Vector of initial values must have length "+ strofreal(parlen)
		start_param = .
	}
	if(rows(start_param) < 1 || start_param == .) {
		
		if(!quiet){
			"Run estimation of three-part nested ordered probit (NOP) model with endogenous switching"
			"Computing starting values:"
		}
		class CNOPModel scalar initial_model 
		initial_model = estimateNOP(y, x, zp, zn, infcat, quiet, ., ., ., ., lambda, maxiter, ptol, vtol, nrtol)
		if(!quiet) {
			"Estimation of NOP with exogenous switching successful"
		}
		start_param = initial_model.params'
		
		if(!quiet) { 
			"Computing starting values for correlation coefficients:"
		}
		latn = 19 * 2 + 1
		ros = rangen(-0.95, 0.95, latn)
		lik = J(latn, latn, 0)

		for(i=1; i<=latn; i++){
			for(j=1; j<=latn; j++) {
				lik[i,j] = sum(MLnopc((start_param, ros[i], ros[j])', x, zp, zn, q, ncat, infcat_index))
			}
		}
		tmp = 0
		maxindex(rowmax(lik),1,i,tmp)
		maxindex(colmax(lik),1,j,tmp)
		ron = ros[j[1,1]]
		rop = ros[i[1,1]]
		start_param = (start_param, rop, ron)'
		if(!quiet) { 
			"rho(+) = " + strofreal(rop)
			"rho(-) = " + strofreal(ron)
		}
	}
	
	if(!quiet) {
		"Estimating NOP with endogenous switching"
	}
	
	// initially transform arguments to avoid inequality constraints
	_nopc_params(start_param, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = ., rop = ., ron = .)
	coded_param = b \ codeIncreasingSequence(a) \ gp \ codeIncreasingSequence(mup) \ gn \ codeIncreasingSequence(mun) \ logit((rop + 1) * 0.5) \ logit((ron + 1) * 0.5)
	/* if coded params contain mistake, replace them by total zero */
	if (max(coded_param :==.)  > 0){
		coded_param = J(rows(coded_param), cols(coded_param), 0)
		if(!quiet){
			"Initial parameters contain error. Replaced with zero"
		}
	}
	
	estimation_successful = 0
	for (estimation_attempt = 0; estimation_attempt <= 7; estimation_attempt++) {
		if(!quiet) {
			"Attempt " + strofreal(estimation_attempt)
		}
		if (estimation_attempt == 0) {
			initial_coded_param = coded_param
			opt_method = "nr"
		}
		if (estimation_attempt == 1) {
			initial_coded_param = coded_param
			opt_method = "bhhh"
		}
		if (estimation_attempt == 2) {
			initial_coded_param = coded_param
			opt_method = "dfp"
		}
		if (estimation_attempt == 3) {
			initial_coded_param = coded_param
			opt_method = "bfgs"
		}
		if (estimation_attempt == 4) {
			initial_coded_param = coded_param * 0
			opt_method = "nr"
		}
		if (estimation_attempt == 5) {
			initial_coded_param = coded_param * 0
			opt_method = "bhhh"
		}
		if (estimation_attempt == 6) {
			initial_coded_param = coded_param * 0
			opt_method = "dfp"
		}
		if (estimation_attempt == 7) {
			initial_coded_param = coded_param * 0
			opt_method = "bfgs"
		}
		singularHmethod = "hybrid" 
		
		
		S = optimize_init()
		if(quiet){
			optimize_init_tracelevel(S , "none")
			optimize_init_verbose(S, 0)
		}
		optimize_init_argument(S, 1, x)
		optimize_init_argument(S, 2, zp)
		optimize_init_argument(S, 3, zn)
		optimize_init_argument(S, 4, q)
		optimize_init_argument(S, 5, ncat)
		optimize_init_argument(S, 6, infcat_index)
		optimize_init_argument(S, 7, 1)
		/* optimize_init_argument(S, 8, lambda)  total L2 regularization term */
		optimize_init_evaluator(S, &_nopc_optim())
		optimize_init_evaluatortype(S, "gf0")
		optimize_init_conv_maxiter(S, maxiter)
		optimize_init_params(S, initial_coded_param')
		optimize_init_conv_ptol(S, ptol)
		optimize_init_conv_vtol(S, vtol)
		optimize_init_conv_nrtol(S, nrtol)
		optimize_init_singularHmethod(S, singularHmethod)
		optimize_init_conv_warning(S, "off") // show that convergence not achieved
		optimize_init_technique(S, opt_method) 
		if (cols(who) > 0 && who != .) {
			optimize_init_cluster(S, who) // added!!!
		}
		errorcode = _optimize(S)
		retCode	= optimize_result_errortext(S)
		convg	= optimize_result_converged(S)
		params = optimize_result_params(S)'
		iterations = optimize_result_iterations(S)
		if (errorcode == 0) {
			if (convg == 1) {
				estimation_successful = 1
				break; // if optimization is successful, then exit
			} else {
				if(!quiet) {
					"Convergence not achieved"
				}
			}
		} else {
			if(!quiet) {
				"NOPC final optimization encountered error code " + strofreal(errorcode) + ": " + retCode
			}
		}
	}
	if (estimation_successful == 0) {
		"Sorry, but despite all attempts, estimation of NOPC parameters did not converge. Turn verbosity on to see the details."
		"Perhaps, there are too few data for such a complex model."
		"Error code is " + strofreal(errorcode) + ": " + retCode
		"Convergence status is " + strofreal(convg)
	}
	
	// After optimization, transform argument back
	
	_nopc_params(params, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = ., rop = ., ron = .)
	params = b \ decodeIncreasingSequence(a) \ gp \ decodeIncreasingSequence(mup) \ gn \ decodeIncreasingSequence(mun) \ invlogit(rop) * 2 - 1 \ invlogit(ron) * 2 - 1
	
	S2 = optimize_init()
	optimize_init_evaluator(S2, &_nopc_optim())
	optimize_init_argument(S2, 1, x)
	optimize_init_argument(S2, 2, zp)
	optimize_init_argument(S2, 3, zn)
	optimize_init_argument(S2, 4, q)
	optimize_init_argument(S2, 5, ncat)
	optimize_init_argument(S2, 6, infcat_index)
	optimize_init_argument(S2, 7, 0)
	optimize_init_params(S2, params')
	
	optimize_init_evaluatortype(S2, "gf0")
	optimize_init_conv_maxiter(S2, maxiter)
	optimize_init_conv_ptol(S2, ptol)
	optimize_init_conv_vtol(S2, vtol)
	optimize_init_conv_nrtol(S2, nrtol)
	optimize_init_singularHmethod(S2, "hybrid")
	optimize_init_conv_warning(S2, "off")
	optimize_init_technique(S2, "nr") 
	
	errorcode2 = _optimize_evaluate(S2)
	
	if(!quiet){
		"Estimation completed"
		""
	}
	
	if (estimation_successful == 0) {
		/* When estimation is not successful, robust covatiance matrix cannot be calculated, and ordinary covariance matrix is . */
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = covMat
	} else {
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = optimize_result_V_robust(S2)
	}
	// calculate predicted probabilities for each observation
	prob_obs = MLnopc(params, x, zp, zn, q, ncat, infcat_index, 1)
	
	class CNOPModel scalar model 
	model.model_class = "NOPC"
	
	model.n	= n
	model.k	= kx
	model.ncat	= ncat
	model.ncatn = ncatn
	model.ncatp = ncatp
	model.infcat = infcat_index
	model.allcat = allcat
	
	model = describeModel(model, params, covMat, covMat_rob, maxLik, n, q, prob_obs)
	
	model.retCode = retCode
	model.error_code = errorcode
	model.etime = clock(c("current_time"),"hms") - starttime
	model.converged = convg
	model.iterations = iterations
	
	return(model)
}


//
class CNOPModel scalar estimateCNOP(y, x, zp, zn, infcat, |quiet, startvalues, robust, who, corresp, lambda, maxiter, ptol, vtol, nrtol) {
	if (rows(quiet) < 1) {
		quiet = 0
	}
	if (args() < 11 || rows(lambda) < 1 || lambda ==.) {
		external CNOP_GLOBAL_CONST_LAMBDA
		lambda = CNOP_GLOBAL_CONST_LAMBDA
	}
	if (args() < 12 || rows(maxiter) < 1 || maxiter ==.) {
		external CNOP_GLOBAL_CONST_MAXITER
		maxiter = CNOP_GLOBAL_CONST_MAXITER
	}
	if (args() < 13 || rows(ptol) < 1 || ptol ==.) {
		external CNOP_GLOBAL_CONST_PTOL
		ptol = CNOP_GLOBAL_CONST_PTOL
	}
	if (args() < 14 || rows(vtol) < 1 || vtol ==.) {
		external CNOP_GLOBAL_CONST_VTOL
		vtol = CNOP_GLOBAL_CONST_VTOL
	}
	if (args() < 15 || rows(nrtol) < 1 || nrtol ==.) {
		external CNOP_GLOBAL_CONST_NRTOL
		nrtol = CNOP_GLOBAL_CONST_NRTOL
	}
	
	starttime = clock(c("current_time"),"hms")
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	allcat = uniqrows(y)
	ncat = rows(allcat)
	ncatp = sum(allcat :> infcat)
	ncatn = sum(allcat :< infcat)
	infcat_index = selectindex(allcat :== infcat)
	
	// incopmlete!!! fill all parameters as in ZIOP
	
	// compute categories
	q = J(n, ncat, 0)
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== allcat[i])
	}
	q0 = (y:==infcat)
	qp = (y:>infcat)
	qn = (y:<infcat)
	q3 = q0 , qp, qn
	
	// initial values
	start_param	= startvalues //  avoid overwriting external argument
	
	parlen =  (kx+2+kzp+ncatp+kzn+ncatn)
	if (rows(start_param) < parlen && rows(start_param) > 0 && start_param != .) {
		"Vector of initial values must have length "+ strofreal(parlen)
		start_param = .
	}
	if(rows(start_param) < 1 || start_param == .) {
		if (!quiet){
			"Computing starting values:"
			"Independent regime equation"
		}
		//start_mu1	= invnormal(runningsum(mean(q3))[1::2])';
		//start_b1		= invsym(x'*x)*x'*(qp-qn);
		//start_param1	= start_b1 \ start_mu1
		params1 = coeffOP(x, q3, 3, quiet, ., lambda, maxiter, ptol, vtol, nrtol)
		
		fltr = (y:>=infcat)
		if (!quiet){
			"Independent outcome equation for y>=0"
		}
		y2 = select(y, fltr)
		q2 = select(q, fltr)
		q2 = select(q2, (J(1,ncatn,0),J(1,ncatp+1,1))) 
		zp2 = select(zp, fltr)
		
		//start_mu1	= invnormal(runningsum(mean(q2))[1::ncatp])';
		//start_b1		= invsym(zp2'*zp2)*zp2'*y2;
		//start_param1	= start_b1 \ start_mu1 
		
		params2 = coeffOP(zp2, q2, ncatp+1,quiet, ., lambda, maxiter, ptol, vtol, nrtol)
		if (!quiet) {
			"Independent outcome equation for y<=0"
		}
		fltr = (y:<=infcat)
		y2 = select(y, fltr)
		q2 = select(q, fltr)
		q2 = select(q2, (J(1,ncatn+1,1),J(1,ncatp,0))) 
		zn2 = select(zn, fltr)
		
		//start_mu1	= invnormal(runningsum(mean(q2))[1::ncatn])';
		//start_b1		= invsym(zn2'*zn2)*zn2'*y2;
		//start_param1	= start_b1 \ start_mu1
		
		params3 = coeffOP(zn2, q2, ncatn+1,quiet, ., lambda, maxiter, ptol, vtol, nrtol)
		
		start_param = params1 \ params2 \ params3
	}
	if(!quiet){
		"ZIOP-3 with exogenous switching"
	}
	
	// initially transform arguments to avoid inequality constraints
	_cnop_params(start_param, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = .)
	coded_param = b \ codeIncreasingSequence(a) \ gp \ codeIncreasingSequence(mup) \ gn \ codeIncreasingSequence(mun)
	/* if coded params contain mistake, replace them by total zero */
	if (max(coded_param :==.)  > 0){
		coded_param = J(rows(coded_param), cols(coded_param), 0)
		if(!quiet){
			"Initial parameters contain error. Replaced with zero"
		}
	}
	
	// make a few attempts to estimate parameters;
	estimation_successful = 0
	for (estimation_attempt = 0; estimation_attempt <= 7; estimation_attempt++) {
		if(!quiet) {
			"Attempt " + strofreal(estimation_attempt)
		}
		if (estimation_attempt == 0) {
			initial_coded_param = coded_param
			opt_method = "nr"
		}
		if (estimation_attempt == 1) {
			initial_coded_param = coded_param
			opt_method = "bhhh"
		}
		if (estimation_attempt == 2) {
			initial_coded_param = coded_param
			opt_method = "dfp"
		}
		if (estimation_attempt == 3) {
			initial_coded_param = coded_param
			opt_method = "bfgs"
		}
		if (estimation_attempt == 4) {
			initial_coded_param = coded_param * 0
			opt_method = "nr"
		}
		if (estimation_attempt == 5) {
			initial_coded_param = coded_param * 0
			opt_method = "bhhh"
		}
		if (estimation_attempt == 6) {
			initial_coded_param = coded_param * 0
			opt_method = "dfp"
		}
		if (estimation_attempt == 7) {
			initial_coded_param = coded_param * 0
			opt_method = "bfgs"
		}
		//singularHmethod = "m-marquardt" 
		singularHmethod = "hybrid" // this one always leads to better (although slower) convergence
		
		
		
		S = optimize_init()
		if(quiet){
			optimize_init_tracelevel(S , "none")
			optimize_init_verbose(S, 0)
		}
		
		optimize_init_argument(S, 1, x)
		optimize_init_argument(S, 2, zp)
		optimize_init_argument(S, 3, zn)
		optimize_init_argument(S, 4, q)
		optimize_init_argument(S, 5, ncat)
		optimize_init_argument(S, 6, infcat_index)
		optimize_init_evaluator(S, &_cnop_optim())
		optimize_init_evaluatortype(S, "gf0") // not concave if use GF1
		optimize_init_conv_maxiter(S, maxiter)
		optimize_init_conv_ptol(S, ptol)
		optimize_init_conv_vtol(S, vtol)
		optimize_init_conv_nrtol(S, nrtol)
		optimize_init_singularHmethod(S, singularHmethod)
		optimize_init_conv_warning(S, "off") // show that convergence not achieved
		optimize_init_technique(S, opt_method) 
		
		if (cols(who) > 0 && who != .) {
			optimize_init_cluster(S, who) // added!!!
		}
		//optimize_init_tracelevel(S, "step")
		
		optimize_init_argument(S, 7, 1) // flag that params are coded to avoid inequality constraints
		optimize_init_params(S, initial_coded_param')
		errorcode = _optimize(S)
		retCode	= optimize_result_errortext(S)
		convg	= optimize_result_converged(S)
		params = optimize_result_params(S)'
		iterations = optimize_result_iterations(S)
		
		if (errorcode == 0) {
			if (convg == 1) {
				estimation_successful = 1
				break; // if optimization is successful, then exit
			} else {
				if(!quiet) {
					"Convergence not achieved"
				}
			}
		} else {
			if(!quiet) {
				"ZIOP-3 final optimization encountered error code " + strofreal(errorcode) + ": " + retCode
			}
		}
	}
	if (estimation_successful == 0) {
		"Sorry, but despite all attempts, estimation of CNOP parameters did not converge. Turn verbosity on to see the details."
		"Perhaps, there are too few data for such a complex model."
		"Error code is " + strofreal(errorcode) + ": " + retCode
		"Convergence status is " + strofreal(convg)
	}
	
	// After optimization, transform argument back
	_cnop_params(params, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = .)
	params = b \ decodeIncreasingSequence(a) \ gp \ decodeIncreasingSequence(mup) \ gn \ decodeIncreasingSequence(mun)
	
	S2 = optimize_init()
	optimize_init_evaluator(S2, &_cnop_optim()) // because of unknown bug I have to repeat this command, otherwise start params are not reset
	optimize_init_argument(S2, 1, x)
	optimize_init_argument(S2, 2, zp)
	optimize_init_argument(S2, 3, zn)
	optimize_init_argument(S2, 4, q)
	optimize_init_argument(S2, 5, ncat)
	optimize_init_argument(S2, 6, infcat_index)
	
	optimize_init_evaluatortype(S2, "gf0")
	optimize_init_conv_maxiter(S2, maxiter)
	optimize_init_conv_ptol(S2, ptol)
	optimize_init_conv_vtol(S2, vtol)
	optimize_init_conv_nrtol(S2, nrtol)
	optimize_init_singularHmethod(S2, "hybrid")
	optimize_init_conv_warning(S2, "off") // show that convergence not achieved
	optimize_init_technique(S2, "nr") 
	
	
	optimize_init_argument(S2, 7, 0)
	optimize_init_params(S2, params')
	//" DDCheck: evaluate start "
	errorcode2 = _optimize_evaluate(S2)
	//" DDCheck: evaluate end "
	
	if(!quiet){
		"Estimation completed"
		""
	}
	
	if (estimation_successful == 0) {
		/* When estimation is not successful, robust covatiance matrix cannot be calculated, and ordinary covariance matrix is . */
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = covMat
	} else {
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = optimize_result_V_robust(S2)
	}
	// calculate predicted probabilities for each observation
	prob_obs = MLcnop(params, x, zp, zn, q, ncat, infcat_index, 1)
	
	class CNOPModel scalar model 
	model.model_class = "CNOP"
	
	model.n	= n
	model.k	= kx
	model.ncat	= ncat
	model.ncatn = ncatn
	model.ncatp = ncatp
	model.infcat = infcat_index
	model.allcat = allcat
	
	model = describeModel(model, params, covMat, covMat_rob, maxLik, n, q, prob_obs)
	
	model.retCode = retCode
	model.error_code = errorcode
	model.etime = clock(c("current_time"),"hms") - starttime
	model.converged = convg
	model.iterations = iterations
	
	return(model)
}


// workfunction that returns CNOPC model (as an instance of CNOPModel class)
class CNOPModel scalar estimateCNOPC(y, x, zp, zn, infcat,|quiet, startvalues, robust, who, corresp, lambda, maxiter, ptol, vtol, nrtol){
	if (rows(quiet) < 1) {
		quiet = 0
	}
	if (args() < 11 || rows(lambda) < 1 || lambda ==.) {
		external CNOP_GLOBAL_CONST_LAMBDA
		lambda = CNOP_GLOBAL_CONST_LAMBDA
	}
	if (args() < 12 || rows(maxiter) < 1 || maxiter ==.) {
		external CNOP_GLOBAL_CONST_MAXITER
		maxiter = CNOP_GLOBAL_CONST_MAXITER
	}
	if (args() < 13 || rows(ptol) < 1 || ptol ==.) {
		external CNOP_GLOBAL_CONST_PTOL
		ptol = CNOP_GLOBAL_CONST_PTOL
	}
	if (args() < 14 || rows(vtol) < 1 || vtol ==.) {
		external CNOP_GLOBAL_CONST_VTOL
		vtol = CNOP_GLOBAL_CONST_VTOL
	}
	if (args() < 15 || rows(nrtol) < 1 || nrtol ==.) {
		external CNOP_GLOBAL_CONST_NRTOL
		nrtol = CNOP_GLOBAL_CONST_NRTOL
	}
	starttime = clock(c("current_time"),"hms")
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	allcat = uniqrows(y)
	ncat = rows(allcat)
	ncatp = sum(allcat:>infcat)
	ncatn = sum(allcat:< infcat)
	infcat_index = selectindex(allcat :== infcat)
	// incopmlete!!! fill all parameters as in ZIOP
	
	// compute categories
	q = J(n, ncat, 0)
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== allcat[i])
	}
	q0 = (y:==infcat)
	qp = (y:>infcat)
	qn = (y:<infcat)
	q3 = q0 , qp, qn
	
	// initial values
	start_param	= startvalues //  avoid overwriting external argument
	
	parlen =  (kx+2+kzp+ncatp+kzn+ncatn+2)
	if (rows(start_param) < parlen && rows(start_param) > 0 && start_param != .) {
		"Vector of initial values must have length "+ strofreal(parlen)
		start_param = .
	}
	if(rows(start_param) < 1 || start_param == .) {
		
		if(!quiet) {
			"Run estimation of three-part zero-inflated ordered probit (ZIOP-3) model with endogenous switching"
		}
		class CNOPModel scalar initial_model 
		initial_model = estimateCNOP(y, x, zp, zn, infcat, quiet, ., ., ., ., lambda, maxiter, ptol, vtol, nrtol)
		if(!quiet) {
			"Estimation of ZIOP-3 with exogenous switching successful"
		}
		start_param = initial_model.params'
		
		if(!quiet) { 
			"Computing starting values for correlation coefficients:"
		}
		latn = 19 * 2 + 1
		ros = rangen(-0.95,0.95,latn)
		lik = J(latn, latn, 0)
		for(i=1; i<=latn; i++){
			for(j=1; j<=latn; j++) {
				lik[i,j] = sum(MLcnopc((start_param, ros[i], ros[j])', x, zp, zn, q, ncat, infcat_index))
			}
		}
		tmp = 0
		maxindex(rowmax(lik),1,i,tmp)
		maxindex(colmax(lik),1,j,tmp)
		rop = ros[i[1,1]]
		ron = ros[j[1,1]]
		start_param = (start_param, rop, ron)'
		if(!quiet) { 
			"rho(+) = " + strofreal(rop)
			"rho(-) = " + strofreal(ron)
		}
	}
	
	if(!quiet) {
		"Estimating ZIOP-3 with endogenous switching"
	}
	
	// initially transform arguments to avoid inequality constraints
	_cnopc_params(start_param, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = ., rop = ., ron = .)
	coded_param = b \ codeIncreasingSequence(a) \ gp \ codeIncreasingSequence(mup) \ gn \ codeIncreasingSequence(mun) \ logit((rop + 1) * 0.5) \ logit((ron + 1) * 0.5)
	/* if coded params contain mistake, replace them by total zero */
	/* INCOMPLETE: need to do it before seeking ros */
	if (max(coded_param :==.)  > 0){
		coded_param = J(rows(coded_param), cols(coded_param), 0)
		if(!quiet){
			"Initial parameters contain error. Replaced with zero"
		}
	}
	
	
	estimation_successful = 0
	for (estimation_attempt = 0; estimation_attempt <= 7; estimation_attempt++) {
		if(!quiet) {
			"Attempt " + strofreal(estimation_attempt)
		}
		if (estimation_attempt == 0) {
			initial_coded_param = coded_param
			opt_method = "nr"
		}
		if (estimation_attempt == 1) {
			initial_coded_param = coded_param
			opt_method = "bhhh"
		}
		if (estimation_attempt == 2) {
			initial_coded_param = coded_param
			opt_method = "dfp"
		}
		if (estimation_attempt == 3) {
			initial_coded_param = coded_param
			opt_method = "bfgs"
		}
		if (estimation_attempt == 4) {
			initial_coded_param = coded_param * 0
			opt_method = "nr"
		}
		if (estimation_attempt == 5) {
			initial_coded_param = coded_param * 0
			opt_method = "bhhh"
		}
		if (estimation_attempt == 6) {
			initial_coded_param = coded_param * 0
			opt_method = "dfp"
		}
		if (estimation_attempt == 7) {
			initial_coded_param = coded_param * 0
			opt_method = "bfgs"
		}
		//singularHmethod = "m-marquardt" 
		singularHmethod = "hybrid" // this one always leads to better (although slower) convergence
		
		
		S = optimize_init()
		if(quiet){
			optimize_init_tracelevel(S , "none")
			optimize_init_verbose(S, 0)
		}
		optimize_init_argument(S, 1, x)
		optimize_init_argument(S, 2, zp)
		optimize_init_argument(S, 3, zn)
		optimize_init_argument(S, 4, q)
		optimize_init_argument(S, 5, ncat)
		optimize_init_argument(S, 6, infcat_index)
		optimize_init_argument(S, 7, 1)
		optimize_init_argument(S, 8, lambda) /* total L2 regularization term */
		optimize_init_evaluator(S, &_cnopc_optim())
		optimize_init_evaluatortype(S, "gf0")
		optimize_init_conv_maxiter(S, maxiter)
		optimize_init_params(S, initial_coded_param')
		
		optimize_init_singularHmethod(S, singularHmethod)
		optimize_init_conv_warning(S, "off") // show that convergence not achieved
		optimize_init_technique(S, opt_method) 
		
		if (cols(who) > 0 && who != .) {
			optimize_init_cluster(S, who) // added!!!
		}

		optimize_init_conv_ptol(S, ptol)
		optimize_init_conv_vtol(S, vtol)
		optimize_init_conv_nrtol(S, nrtol)
		errorcode = _optimize(S)
		
		
		retCode	= optimize_result_errortext(S)
		convg	= optimize_result_converged(S)
		params = optimize_result_params(S)'
		iterations = optimize_result_iterations(S)
		
		if (errorcode == 0) {
			if (convg == 1) {
				estimation_successful = 1
				break; // if optimization is successful, then exit
			} else {
				if(!quiet) {
					"Convergence not achieved"
				}
			}
		} else {
			if(!quiet) {
				"ZIOP-3  final optimization encountered error code " + strofreal(errorcode) + ": " + retCode
			}
		}
	}
	if (estimation_successful == 0) {
		"Sorry, but despite all attempts, estimation of CNOPC parameters did not converge. Turn verbosity on to see the details."
		"Perhaps, there are too few data for such a complex model."
		"Error code is " + strofreal(errorcode) + ": " + retCode
		"Convergence status is " + strofreal(convg)
	}

	// After optimization, transform argument back
	_cnopc_params(params, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = ., rop = ., ron = .)
	params = b \ decodeIncreasingSequence(a) \ gp \ decodeIncreasingSequence(mup) \ gn \ decodeIncreasingSequence(mun) \ invlogit(rop) * 2 - 1 \ invlogit(ron) * 2 - 1
	
	
	S2 = optimize_init()
	optimize_init_evaluator(S2, &_cnopc_optim()) // because of unknown bug I have to repeat this command, otherwise start params are not reset
	
	optimize_init_argument(S2, 1, x)
	optimize_init_argument(S2, 2, zp)
	optimize_init_argument(S2, 3, zn)
	optimize_init_argument(S2, 4, q)
	optimize_init_argument(S2, 5, ncat)
	optimize_init_argument(S2, 6, infcat_index)
	optimize_init_argument(S2, 7, 0)
	optimize_init_argument(S2, 8, 0)
	optimize_init_evaluatortype(S2, "gf0")
	optimize_init_conv_maxiter(S2, maxiter)
	optimize_init_conv_ptol(S2, ptol)
	optimize_init_conv_vtol(S2, vtol)
	optimize_init_conv_nrtol(S2, nrtol)
	optimize_init_singularHmethod(S2, "hybrid")
	optimize_init_conv_warning(S2, "off") // show that convergence not achieved
	optimize_init_technique(S2, "nr") 
	optimize_init_params(S2, params')
	//" DDCheck: evaluate start "
	errorcode = _optimize_evaluate(S2)
	//" DDCheck: evaluate end "
	
	// estimate main miopr stage
	if(!quiet) {
		"Estimation completed"
		""
	}
	// incomplete: bad convergence even at good starting values. CHECK CONVEXITY
	
	if (estimation_successful == 0) {
		/* When estimation is not successful, robust covatiance matrix cannot be calculated, and ordinary covariance matrix is . */
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = covMat
	} else {
		maxLik	= optimize_result_value(S2)
		grad 	= optimize_result_gradient(S2)
		covMat	= optimize_result_V(S2)
		covMat_rob = optimize_result_V_robust(S2)
	}
	// calculate predicted probabilities for each observation
	prob_obs = MLcnopc(params, x, zp, zn, q, ncat, infcat_index, 1)
	
	class CNOPModel scalar model 
	model.model_class = "CNOPC"
	
	model.n	= n
	model.k	= kx
	model.ncat	= ncat
	model.ncatn = ncatn
	model.ncatp = ncatp
	model.infcat = infcat_index
	model.allcat = allcat
	
	model = describeModel(model, params, covMat, covMat_rob, maxLik, n, q, prob_obs)

	model.retCode = retCode
	model.error_code = errorcode
	model.etime = clock(c("current_time"),"hms") - starttime
	model.converged = convg
	model.iterations = iterations
	
	return(model)
}



// for each element of candidates, return its position in subset or 0
function positionsInList(candidates, subset) {
	ans = J(1, cols(candidates), 0)
	for(i = 1; i <= cols(subset); i++) {
		filter = candidates :== subset[i]
		if(sum(filter) > 0){
			ans[selectindex(filter)] = i
		}
	}
	return(ans)
}

// passes CNOP and CNOP(c) specification from Stata to Mata and back
function processCNOP(yxnames, zpnames, znnames, infcat, correlated, touse, robust, cluster, initial, nolog) {
	xytokens = tokens(yxnames)
	yname = xytokens[1]
	xnames = invtokens(xytokens[,2::cols(xytokens)])
	if(strlen(zpnames) == 0) {
		if(strlen(znnames) == 0) {
			znnames = xnames
			zpnames = xnames
		} else {
			zpnames = znnames
		}
	} else if (strlen(znnames) == 0){
		znnames = zpnames
	}
	
	allvars = uniqrows( (tokens(xnames),tokens(zpnames),tokens(znnames))' )'
	corresp = (positionsInList(allvars, tokens(xnames)) \ positionsInList(allvars, tokens(zpnames)) \ positionsInList(allvars, tokens(znnames)) )
	st_view(y  = ., ., yname, touse)
	st_view(x  = ., ., xnames, touse)
	st_view(zp = ., ., zpnames, touse)
	st_view(zn = ., ., znnames, touse)
	who = .
	if (cluster != "") {
		st_view(who = ., ., cluster, touse)
		robust = 1
	}
	
	if (initial != "") {
		initial = strtoreal(tokens(initial))'
		if (sum(initial :== .) > 0) {
			"Incorrect initial values! Expected a numeric sequence delimited with whitespace."
			"Default initial values will be used."
			initial = .
		}
	} else {
		initial = .
	}
	
	if ( (sum(uniqrows(y) :< infcat) < 1) | (sum(uniqrows(y) :> infcat) < 1) ) {
		errprintf("The dependent variable takes on less than three discrete values. ")
		errprintf("The ZIOP-3 model is designed for a dependent variable with at least three outcome choices.")
		exit(3498)
	}
	
	class CNOPModel scalar model 
	if (correlated) {
		switching_type = "endogenous"
		model = estimateCNOPC(y, x, zp, zn, infcat, nolog, initial, robust, who)
	} else {
		switching_type = "exogenous"
		model = estimateCNOP(y, x, zp, zn, infcat, nolog, initial, robust, who)
	}
	
	model.yname = yname
	model.xnames = xnames
	model.zpnames = zpnames
	model.znnames = znnames
	
	model.XZnames 	= allvars
	model.corresp	= corresp
	st_view(xz = ., ., invtokens(allvars))
	model.XZmeans = mean(xz)
	
	model.eqnames = J(1, cols(tokens(xnames)) + 2, "Regime equation"),  J(1, cols(tokens(zpnames)) + model.ncatp, "Outcome equation (+)"),  J(1, cols(tokens(znnames)) + model.ncatn, "Outcome equation (-)")
	model.parnames = tokens(xnames), "/cut1", "/cut2", tokens(zpnames),  "/cut" :+ strofreal(1..model.ncatp), tokens(znnames),  "/cut" :+ strofreal(1..model.ncatn)
	
	if (correlated) {
		model.eqnames = model.eqnames, J(1, 2, "Correlation coefficients")
		model.parnames = model.parnames, "rho(+)", "rho(-)"
	}
	st_matrix("b", model.params')
	if (robust == 1) {
		st_matrix("V", model.V_rob)
	} else {
		st_matrix("V", model.V)
	}
	model.robust = robust
	stripes = model.eqnames' , model.parnames'
	st_matrixcolstripe("b", stripes)
	st_matrixcolstripe("V", stripes)
	st_matrixrowstripe("V", stripes)
	st_local("depvar", yname)
	st_local("N", strofreal(model.n))
	st_local("ll", strofreal(model.logLik))
	st_matrix("ll_obs", model.ll_obs)
	
	// describe the model
	"Three-part zero-inflated ordered probit model with " + switching_type + " switching"
	"Number of observations = " + strofreal(model.n)
	"Log likelihood = " + strofreal(model.logLik)
	"AIC            = " + strofreal(model.AIC)
	"BIC            = " + strofreal(model.BIC)
	
	return(model)
}

function processNOP(yxnames, zpnames, znnames, infcat, correlated, touse, robust, cluster, initial, nolog) {
	xytokens = tokens(yxnames)
	yname = xytokens[1]
	xnames = invtokens(xytokens[,2::cols(xytokens)])
	if(strlen(zpnames) == 0) {
		if(strlen(znnames) == 0) {
			znnames = xnames
			zpnames = xnames
		} else {
			zpnames = znnames
		}
	} else if (strlen(znnames) == 0){
		znnames = zpnames
	}
	
	allvars = uniqrows( (tokens(xnames),tokens(zpnames),tokens(znnames))' )'
	corresp = (positionsInList(allvars, tokens(xnames)) \ positionsInList(allvars, tokens(zpnames)) \ positionsInList(allvars, tokens(znnames)) )
	st_view(y  = ., ., yname, touse)
	st_view(x  = ., ., xnames, touse)
	st_view(zp = ., ., zpnames, touse)
	st_view(zn = ., ., znnames, touse)
	who = .
	if (cluster != "") {
		st_view(who = ., ., cluster, touse)
		robust = 1
	}
	if (initial != "") {
		initial = strtoreal(tokens(initial))'
		if (sum(initial :== .) > 0) {
			"Incorrect initial values! Expected a numeric sequence delimited with whitespace."
			"Default initial values will be used."
			initial = .
		}
	} else {
		initial = .
	}
	
	if ( (sum(uniqrows(y) :< infcat) < 2) | (sum(uniqrows(y) :> infcat) < 2) ) {
		errprintf("The dependent variable takes on less than five discrete values. ")
		errprintf("The NOP model is designed for a dependent variable with at least five outcome choices.")
		errprintf("With only three or two outcome choices the NOP model reduces to the conventional ordered probit model. ")
		errprintf("Use the oprobit command. ")
		exit(3498)
	}
	
	class CNOPModel scalar model 
	if (correlated) {
		switching_type = "endogenous"
		model = estimateNOPC(y, x, zp, zn, infcat, nolog, initial, robust, who)
	} else {
		switching_type = "exogenous"
		model = estimateNOP(y, x, zp, zn, infcat, nolog, initial, robust, who)
	}
	
	model.yname = yname
	model.xnames = xnames
	model.zpnames = zpnames
	model.znnames = znnames
	
	model.XZnames 	= allvars
	model.corresp	= corresp
	st_view(xz = ., ., invtokens(allvars))
	model.XZmeans = mean(xz)
	
	model.eqnames = J(1, cols(tokens(xnames))+2, "Regime equation"),  J(1, cols(tokens(zpnames)) + model.ncatn - 1, "Outcome equation (+)"),  J(1, cols(tokens(znnames)) + model.ncatn - 1, "Outcome equation (-)")
	model.parnames = tokens(xnames), "/cut1", "/cut2", tokens(zpnames),  "/cut" :+ strofreal(1..(model.ncatp-1)), tokens(znnames),  "/cut" :+ strofreal(1..(model.ncatn-1))
	
	if (correlated) {
		model.eqnames = model.eqnames, J(1, 2, "Correlation coefficients")
		model.parnames = model.parnames, "rho(+)", "rho(-)"
	}
	st_matrix("b", model.params')
	if (robust == 1) {
		st_matrix("V", model.V_rob)
	} else {
		st_matrix("V", model.V)
	}
	model.robust = robust
	stripes = model.eqnames' , model.parnames'
	st_matrixcolstripe("b", stripes)
	st_matrixcolstripe("V", stripes)
	st_matrixrowstripe("V", stripes)
	st_local("depvar", yname)
	st_local("N", strofreal(model.n))
	st_local("ll", strofreal(model.logLik))
	st_matrix("ll_obs", model.ll_obs)
	
	// describe the model
	"Three-part nested ordered probit model with " + switching_type + " switching"
	"Number of observations = " + strofreal(model.n)
	"Log likelihood = " + strofreal(model.logLik)
	"AIC            = " + strofreal(model.AIC)
	"BIC            = " + strofreal(model.BIC)
	
	return(model)
}


// passes ZIOP specification from Stata to Mata and back
function processMIOPR(yxnames, znames, infcat, correlated, touse, robust, cluster, initial, nolog) {
	xytokens = tokens(yxnames)
	yname = xytokens[1]
	xnames = invtokens(xytokens[,2::cols(xytokens)])
	if (strlen(znames) == 0){
		znames = xnames
	}
	allvars = uniqrows( (tokens(xnames),tokens(znames))' )'
	corresp = (positionsInList(allvars, tokens(xnames)) \ positionsInList(allvars, tokens(znames)))
	st_view(y  = ., ., yname, touse)
	st_view(x  = ., ., xnames, touse)
	st_view(z = ., ., znames, touse)
	who = .
	if (cluster != "") {
		st_view(who = ., ., cluster, touse)
		robust = 1
	}
	if (initial != "") {
		initial = strtoreal(tokens(initial))'
		if (sum(initial :== .) > 0) {
			"Incorrect initial values! Expected a numeric sequence delimited with whitespace."
			"Default initial values will be used."
			initial = .
		}
	} else {
		initial = .
	}
	
	class CNOPModel scalar model
	if (correlated) {
		switching_type = "endogenous"
		model = estimateMIOPRC(y, x, z, infcat, nolog, initial, robust, who)
	} else {
		switching_type = "exogenous"
		model = estimateMIOPR(y, x, z, infcat, nolog, initial, robust, who)
	}
	model.XZnames 	= allvars
	model.yname = yname
	model.xnames = xnames
	model.znames = znames
	model.corresp	= corresp
	st_view(xz = ., ., invtokens(allvars))
	model.XZmeans = mean(xz)
	model.eqnames = J(1, cols(tokens(xnames)) + 1, "Regime equation"), J(1, cols(tokens(znames)) + rows(model.allcat)-1, "Outcome equation")
	model.parnames = tokens(xnames), "/cut1", tokens(znames), "/cut" :+ strofreal(1..(rows(model.allcat)-1))
	
	if (correlated) {
		model.eqnames = model.eqnames, "Correlation coefficient"
		model.parnames = model.parnames, "rho"
	}
	
	st_matrix("b", model.params')
	if (robust == 1) {
		st_matrix("V", model.V_rob)
	} else {
		st_matrix("V", model.V)
	}
	model.robust = robust
	//kk = rows(model.params)
	//stripes = J(kk, 1, ""), "coef" :+ strofreal(1::kk)
	stripes = model.eqnames' , model.parnames'
	st_matrixcolstripe("b", stripes)
	st_matrixcolstripe("V", stripes)
	st_matrixrowstripe("V", stripes)
	st_local("depvar", yname)
	st_local("N", strofreal(model.n))
	st_local("ll", strofreal(model.logLik))
	st_matrix("ll_obs", model.ll_obs)
	
	// describe the model
	"Two-part zero-inflated ordered probit model with " + switching_type + " switching"
	"Number of observations = " + strofreal(model.n)
	"Log likelihood = " + strofreal(model.logLik)
	"AIC            = " + strofreal(model.AIC)
	"BIC            = " + strofreal(model.BIC)
	
	return(model)
}


function get_colstripes(model_class, loop, allcat, infcat) {
	if (loop == 1) {
		colstripes = "Pr(y=" :+ strofreal(allcat) :+ ")"
	}
	if (loop == 2) {
		if (model_class == "MIOPR" || model_class == "MIOPRC") {
			colstripes = ("Pr(y=0|s=0)" \ "Pr(y=0|s=1)")
		} else {
			colstripes = ("Pr(y=0|s=0)" \ "Pr(y=0|s=-1)" \  "Pr(y=0|s=+1)")
		}
	}
	if (loop == 3) {
		if (model_class == "MIOPR" || model_class == "MIOPRC") {
			colstripes = ("Pr(s=0)" \ "Pr(s=1)")
		} else {
			colstripes = ("Pr(s=-1)" \ "Pr(s=0)" \  "Pr(s=+1)")
		}
	}
	// workaround: stata does not allow colstripes containing dots
	colstripes = subinstr(colstripes, "=.", "=0.")
	colstripes = subinstr(colstripes, "=-.", "=-0.")
	colstripes = subinstr(colstripes, ".", ",")
	return(colstripes)
}

function output_matrix(matrix_name, matrix_value, rowstripes, colstripes){
	st_matrix(matrix_name, matrix_value)
	st_matrixrowstripe(matrix_name, (J(rows(rowstripes), 1, ""), rowstripes))
	st_matrixcolstripe(matrix_name, (J(rows(colstripes), 1, ""), colstripes))
}

function output_mesetp(me, se, rowstripes, colstripes) {
	t = me :/ se
	pval = (1:-normal(abs(t))) :* 2
	output_matrix("r(me)",     me, rowstripes, colstripes)
	output_matrix("r(se)",     se, rowstripes, colstripes)
	output_matrix("r(t)",       t, rowstripes, colstripes)
	output_matrix("r(pval)", pval, rowstripes, colstripes)
}

function update_named_vector(values, names, tokens) {
	atVarnames = tokens[range(1, cols(tokens)-2, 3)]
	atValues = subinstr(tokens[range(3, cols(tokens), 3)], ",", "") 
	atTable = sort( (atVarnames' , atValues'), 1)
	for(i = 1; i <= rows(atTable); i++) {
		index = selectindex(names :== atTable[i,1])
		if (cols(index)) {
			newValue = strtoreal(atTable[i,2])
			if (newValue == .) {
				"'" + atTable[i,1] + " = " + atTable[i,2] + " could not be parsed"
			}
			values[index] = newValue
		} else {
			atTable[i,1] + " was not applied in the last ZIOP/NOP model"
		}
	}
	return(values)
} 

// marginal effects for MIOP(r), CNOP, CNOP(c)

function CNOPmargins(class CNOPModel scalar model, string atVarlist, string dummiesVarlist, zeroes, regime) {
	dummiesVector = positionsInList(model.XZnames, tokens(dummiesVarlist))
	xzbar = model.XZmeans
	atTokens = tokens(atVarlist, " =")
	
	if (length(atTokens) >= 3) {
		xzbar = update_named_vector(xzbar, model.XZnames, atTokens)
	}
	loop = 1 // code of prediction type
	if (zeroes) {
		loop = 2
	} else if (regime) {
		loop = 3
	}
	output_matrix("r(at_all)", xzbar, " ", model.XZnames')
	
	rowstripes = model.XZnames'
	colstripes = get_colstripes(model.model_class, loop, model.allcat, model.infcat)
	
	mese = generalMEwithSE(xzbar, model, dummiesVector, loop)
	kxz = cols(xzbar)
	me = mese[1::kxz,]
	se = mese[(1::kxz) :+ kxz,]
	
	output_mesetp(me, se, rowstripes, colstripes)	
}


function CNOPprobabilities(class CNOPModel scalar model, string atVarlist, zeroes, regime) {
	xz_from = model.XZmeans
	atTokens = tokens(atVarlist, " =")
	
	if (length(atTokens) >= 3) {
		xz_from = update_named_vector(xz_from, model.XZnames, atTokens)
	}
	
	loop = 1 // code of prediction type
	if (zeroes) {
		loop = 2
	} else if (regime) {
		loop = 3
	}
	
	output_matrix("r(at_all)", xz_from, " ", model.XZnames')

	colstripes = get_colstripes(model.model_class, loop, model.allcat, model.infcat)
	rowstripes = " " // rowstripes made invisible
	mese = generalPredictWithSE(xz_from, model, loop)
	me = mese[1,]
	se = mese[2,]
	output_mesetp(me, se, rowstripes, colstripes)
}

function CNOPcontrasts(class CNOPModel scalar model, string atVarlist, string toVarlist, zeroes, regime) {
	xz_from = model.XZmeans
	xz_to = model.XZmeans
	atTokens = tokens(atVarlist, " =")
	toTokens = tokens(toVarlist, " =")
	
	if (length(atTokens) >= 3) {
		xz_from = update_named_vector(xz_from, model.XZnames, atTokens)
	}
	
	if (length(toTokens) >= 3) {
		xz_to = update_named_vector(xz_to, model.XZnames, toTokens)
	}
	
	loop = 1 // code of prediction type
	if (zeroes) {
		loop = 2
	} else if (regime) {
		loop = 3
	}
	
	if (sum((xz_from - xz_to):^2) == 0) {
		"Trying to contrast the same point"
	}
	
	output_matrix("r(between_all)", xz_from \ xz_to, "from" \ "to", model.XZnames')
	
	colstripes = get_colstripes(model.model_class, loop, model.allcat, model.infcat)
	rowstripes = " " // rowstripes made invisible
	
	mese = generalContrastsWithSE(xz_from, xz_to, model, loop)
	me = mese[1,]
	se = mese[2,]
	output_mesetp(me, se, rowstripes, colstripes)
}

// prediction for MIOP(r), CNOP, CNOP(c)
function CNOP_predict(class CNOPModel scalar model, string scalar newVarName, real scalar zeroes, real scalar regime, string scalar touse, scalar output) {
	sp = strpos(newVarName, ",")
	if (sp != 0){
		newVarName = substr(newVarName, 1, sp - 1)
	}
	loop = 1 // code of prediction type
	if (zeroes) {
		loop = 2
	} else if (regime) {
		loop = 3
	}
	label_indices = strofreal(1..model.ncat)
	labels = strofreal(model.allcat')
	values = model.allcat'
	if (model.model_class == "MIOPR" || model.model_class == "MIOPRC") {
		st_view(x  = ., ., model.xnames, touse)
		st_view(z  = ., ., model.znames, touse)
		if (model.model_class == "MIOPR") {
			p 	= MLmiopr(model.params, x, z, q=., model.ncat, model.infcat, loop)
		} else {
			p 	= MLmioprc(model.params, x, z, q=., model.ncat, model.infcat, loop)
		}
		if (loop == 2) {
			label_indices = ("0", "1")
			labels = labels = ("zero 0", "non-zero 0")
			values = (0, 1)
		}
		if (loop == 3) {
			label_indices = ("0", "1")
			labels = ("neutral", "active") :+ " regime"
			values = (0, 1)
		}
	} else if (model.model_class == "CNOP" || model.model_class == "CNOPC") {
		st_view(x  = ., ., model.xnames, touse)
		st_view(zp = ., ., model.zpnames, touse)
		st_view(zn = ., ., model.znnames, touse)
		if (model.model_class == "CNOP") {
			p 	= MLcnop(model.params, x, zp, zn, q=., model.ncat, model.infcat, loop)
		} else {
			p 	= MLcnopc(model.params, x, zp, zn, q=., model.ncat, model.infcat, loop)
		}
		if (loop == 2) {
			label_indices = ("0", "n", "p")
			labels = ("zero 0", "negative 0", "positive 0")
			values = (0, -1, 1)
		}
		if (loop == 3) {
			label_indices = ("n", "0", "p")
			labels = ("negative", "neutral",  "positive") :+ " regime"
			values = (-1, 0, 1)
		}
	} else if (model.model_class == "NOP" || model.model_class == "NOPC") {
		st_view(x  = ., ., model.xnames, touse)
		st_view(zp = ., ., model.zpnames, touse)
		st_view(zn = ., ., model.znnames, touse)
		if (model.model_class == "NOP") {
			p 	= MLnop(model.params, x, zp, zn, q=., model.ncat, model.infcat, loop)
		} else {
			p 	= MLnopc(model.params, x, zp, zn, q=., model.ncat, model.infcat, loop)
		}
		if (loop == 2) {
			label_indices = ("0", "n", "p")
			labels = ("zero 0", "negative 0", "positive 0")
			values = (0, -1, 1)
		}
		if (loop == 3) {
			label_indices = ("n", "0", "p")
			labels = ("negative", "neutral",  "positive") :+ " regime"
			values = (-1, 0, 1)
		}
	}
	label_indices = newVarName + "_" :+ label_indices
	if (output == "mode") {
		if (_st_varindex(newVarName) :== .) {
			tmp = st_addvar("double", newVarName)
		}
		st_view(v = ., ., newVarName)
		prediction = rowsum((p:==rowmax(p)) :* values)
		v[,] = prediction
		st_vlmodify(newVarName, values', labels')
	} else if (output == "mean") {
		if (_st_varindex(newVarName) :== .) {
			tmp = st_addvar("double", newVarName)
		}
		st_view(v = ., ., newVarName)
		prediction = rowsum(p :* values)
		v[,] = prediction
		st_vlmodify(newVarName, values', labels')
	} else if (output == "cum"){
		tmp = st_addvar("double", label_indices[selectindex(_st_varindex(label_indices) :== .)])
		st_view(v = ., ., label_indices)
		for (i = 1; i <= length(labels); ++i) {
			st_varlabel(label_indices[i], labels[i])
		}
		v[,1] = p[,1]
		for (i = 2; i <= cols(p); ++i) {
			v[,i] = v[,i-1] + p[,i]
		}
	} else {
		tmp = st_addvar("double", label_indices[selectindex(_st_varindex(label_indices) :== .)])
		st_view(v = ., ., label_indices)
		for (i = 1; i <= length(labels); ++i) {
			st_varlabel(label_indices[i], labels[i])
		}
		v[,] = p
	}
}
end

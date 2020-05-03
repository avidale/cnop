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

function colmedian(input_matrix) {
	m = cols(input_matrix)
	n = rows(input_matrix)
	result = J(1, m, .)
	for (i = 1; i <= m; i++) {
		sorted = sort(input_matrix[ , i], 1)
		result[1, i] = sorted[(n + 1) / 2]
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
	
	model.df = rows(params)
	model.df_null = cols(q) - 1
	model.chi2 = 2 * (maxLik - model.logLik0)
	model.chi2_pvalue = 1 - chi2(model.df - model.df_null, model.chi2)
	
	model.brier_score = matrix_mse(prob_obs - q)
	model.ranked_probability_score = matrix_mse(running_rowsum(prob_obs) - running_rowsum(q))
	
	values = runningsum(J(1, cols(q), 1))
	prediction = rowsum((prob_obs:==rowmax(prob_obs)) :* values)
	actual = rowsum((q:==rowmax(q)) :* values)
	model.accuracy = mean(prediction :== actual)
	
	model.V	= covMat
	model.V_rob	= covMat_rob
	model.logLik	= maxLik
	model.probabilities = prob_obs
	model.ll_obs = log(rowsum(prob_obs :* q))
	
	return (model)
}


// calculate and remember even more postestimation statistics
class CNOPModel scalar postDescribeModel(class CNOPModel scalar model, robust, allvars, corresp) {
	model.robust = robust
	model.XZnames 	= allvars
	model.corresp	= corresp
	st_view(xz = ., ., invtokens(allvars))
	model.XZmeans = mean(xz)
	model.XZmedians = colmedian(xz)
	return (model)
}



// workfunction that returns OP model (as an instance of CNOPModel class)
class CNOPModel scalar estimateOP(y, x, |quiet, startvalues, xbar, dummies, robust, who) {
	
	starttime = clock(c("current_time"), "hms")
	
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
	optimize_init_argument(S, 4, 0)
	optimize_init_evaluator(S, &_op_optim())
	optimize_init_evaluatortype(S, "gf0") // todo: gf1 raises discontinuity error - why ?
	optimize_init_params(S, (start_param'))
	
	if (cols(who) > 0 && who != .) {
		optimize_init_cluster(S, who)
	}
	errorcode = _optimize(S)
	if (errorcode == 0) {
		params = optimize_result_params(S)'
	} else {
		"dd: OP final optimization encountered error code " + strofreal(errorcode)
		params = optimize_result_params(S)'
	}
	
	// extract optimization results
	if (errorcode == 0) {
		/* When estimation is not successful, robust covatiance matrix cannot be calculated, and ordinary covariance matrix is . */
		maxLik	= optimize_result_value(S)
		grad 	= optimize_result_gradient(S)
		covMat	= optimize_result_V(S)
		covMat_rob = covMat
	} else {
		maxLik	= optimize_result_value(S)
		grad 	= optimize_result_gradient(S)
		covMat	= optimize_result_V(S)
		covMat_rob = optimize_result_V_robust(S)
	}
	
	prob_obs	= MLop(params, x, q, ncat, 1)
	
	// pack results
	class CNOPModel scalar model 
	model.model_class = "OP"
	
	model.n	= n
	model.k	= k
	model.ncat	= ncat
	model.allcat = allcat
	
	model = describeModel(model, params, covMat, covMat_rob, maxLik, n, q, prob_obs)
	
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
		lambda = CNOP_GLOBAL_CONST_LAMBDA()
	}
	if (args() < 12 || rows(maxiter) < 1 || maxiter ==.) {
		maxiter = CNOP_GLOBAL_CONST_MAXITER()
	}
	if (args() < 13 || rows(ptol) < 1 || ptol ==.) {
		ptol = CNOP_GLOBAL_CONST_PTOL()
	}
	if (args() < 14 || rows(vtol) < 1 || vtol ==.) {
		vtol = CNOP_GLOBAL_CONST_VTOL()
	}
	if (args() < 15 || rows(nrtol) < 1 || nrtol ==.) {
		nrtol = CNOP_GLOBAL_CONST_NRTOL()
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
	} else {
		if(!quiet) {
			"Using given start values for parameter estimation"
		}
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
		lambda = CNOP_GLOBAL_CONST_LAMBDA()
	}
	if (args() < 11 || rows(maxiter) < 1 || maxiter ==.) {
		maxiter = CNOP_GLOBAL_CONST_MAXITER()
	}
	if (args() < 12 || rows(ptol) < 1 || ptol ==.) {
		ptol = CNOP_GLOBAL_CONST_PTOL()
	}
	if (args() < 13 || rows(vtol) < 1 || vtol ==.) {
		vtol = CNOP_GLOBAL_CONST_VTOL()
	}
	if (args() < 14 || rows(nrtol) < 1 || nrtol ==.) {
		nrtol = CNOP_GLOBAL_CONST_NRTOL()
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
	}  else {
		if(!quiet) {
			"Using given start values for parameter estimation"
		}
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
		lambda = CNOP_GLOBAL_CONST_LAMBDA()
	}
	if (args() < 11 || rows(maxiter) < 1 || maxiter ==.) {
		maxiter = CNOP_GLOBAL_CONST_MAXITER()
	}
	if (args() < 12 || rows(ptol) < 1 || ptol ==.) {
		ptol = CNOP_GLOBAL_CONST_PTOL()
	}
	if (args() < 13 || rows(vtol) < 1 || vtol ==.) {
		vtol = CNOP_GLOBAL_CONST_VTOL()
	}
	if (args() < 14 || rows(nrtol) < 1 || nrtol ==.) {
		nrtol = CNOP_GLOBAL_CONST_NRTOL()
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
	} else {
		if(!quiet) {
			"Using given start values for parameter estimation"
		}
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
		lambda = CNOP_GLOBAL_CONST_LAMBDA()
	}
	if (args() < 12 || rows(maxiter) < 1 || maxiter ==.) {
		maxiter = CNOP_GLOBAL_CONST_MAXITER()
	}
	if (args() < 13 || rows(ptol) < 1 || ptol ==.) {
		ptol = CNOP_GLOBAL_CONST_PTOL()
	}
	if (args() < 14 || rows(vtol) < 1 || vtol ==.) {
		vtol = CNOP_GLOBAL_CONST_VTOL()
	}
	if (args() < 15 || rows(nrtol) < 1 || nrtol ==.) {
		nrtol = CNOP_GLOBAL_CONST_NRTOL()
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
	} else {
		if(!quiet) {
			"Using given start values for parameter estimation"
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
		lambda = CNOP_GLOBAL_CONST_LAMBDA()
	}
	if (args() < 12 || rows(maxiter) < 1 || maxiter ==.) {
		maxiter = CNOP_GLOBAL_CONST_MAXITER()
	}
	if (args() < 13 || rows(ptol) < 1 || ptol ==.) {
		ptol = CNOP_GLOBAL_CONST_PTOL()
	}
	if (args() < 14 || rows(vtol) < 1 || vtol ==.) {
		vtol = CNOP_GLOBAL_CONST_VTOL()
	}
	if (args() < 15 || rows(nrtol) < 1 || nrtol ==.) {
		nrtol = CNOP_GLOBAL_CONST_NRTOL()
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
	} else {
		if(!quiet) {
			"Using given start values for parameter estimation"
		}
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
		lambda = CNOP_GLOBAL_CONST_LAMBDA()
	}
	if (args() < 12 || rows(maxiter) < 1 || maxiter ==.) {
		maxiter = CNOP_GLOBAL_CONST_MAXITER()
	}
	if (args() < 13 || rows(ptol) < 1 || ptol ==.) {
		ptol = CNOP_GLOBAL_CONST_PTOL()
	}
	if (args() < 14 || rows(vtol) < 1 || vtol ==.) {
		vtol = CNOP_GLOBAL_CONST_VTOL()
	}
	if (args() < 15 || rows(nrtol) < 1 || nrtol ==.) {
		nrtol = CNOP_GLOBAL_CONST_NRTOL()
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
	} else {
		if(!quiet) {
			"Using given start values for parameter estimation"
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
	/* return a vector where in i'th position is 0, if i'th candidate is not in subset, otherwise its position in subset */
	ans = J(1, cols(candidates), 0)
	for(i = 1; i <= cols(subset); i++) {
		filter = candidates :== subset[i]
		if(sum(filter) > 0){
			ans[selectindex(filter)] = i
		}
	}
	return(ans)
}

function vuong_vs_op(class CNOPModel scalar model) {
	st_view(op_y=., ., invtokens(model.yname))
	st_view(op_x=., ., invtokens(model.XZnames))
	class CNOPModel scalar op_model 
	op_model = estimateOP(op_y, op_x, 1)
	ll_diff = model.ll_obs - op_model.ll_obs
	k_1 = rows(model.params)
	k_2 = rows(op_model.params)
	vuong_calc(ll_diff, k_1, k_2)
}

//
function passModelToStata(class CNOPModel scalar model) {
	
	if (model.model_class == "NOPC" | model.model_class == "MIOPRC" | model.model_class == "CNOPC") {
		switching_type = "endogenous"
	} else {
		switching_type = "exogenous"
	}
	
	if (model.model_class == "NOP" | model.model_class == "NOPC") {
		model_suptype = "Nested ordered probit regression"
		model_type = "Three-part nested ordered probit model"
		inflation_line = ""
	} else if (model.model_class == "CNOP" | model.model_class == "CNOPC") {
		model_suptype = "Zero-inflated ordered probit regression"
		model_type = "Three-part zero-inflated ordered probit model"
		inflation_line = "Zero inflation:          three regimes"
	} else if (model.model_class == "MIOPR" | model.model_class == "MIOPRC") {
		model_suptype = "Zero-inflated ordered probit regression"
		model_type = "Two-part zero-inflated ordered probit model"
		inflation_line = "Zero inflation:          two regimes"
	}

	st_matrix("b", model.params')
	if (model.robust == 1) {
		st_matrix("V", model.V_rob)
	} else {
		st_matrix("V", model.V)
	}
	stripes = model.eqnames' , model.parnames'
	st_matrixcolstripe("b", stripes)
	st_matrixcolstripe("V", stripes)
	st_matrixrowstripe("V", stripes)
	st_local("depvar", model.yname)
	st_local("N", strofreal(model.n))
	st_numscalar("ll", model.logLik)
	st_numscalar("k", rows(model.params))
	st_matrix("ll_obs", model.ll_obs)
	st_numscalar("r2_p", model.R2)
	st_numscalar("k_cat", model.ncat)
	st_numscalar("df_m", model.df)
	st_numscalar("ll_0", model.logLik0)
	st_numscalar("chi2", model.chi2)
	st_numscalar("p", model.chi2_pvalue)
	st_numscalar("aic", model.AIC)
	st_numscalar("bic", model.BIC)
	
	// describe the model
	//model_type + " with " + switching_type + " switching"
	displayas("result")
	printf("%s\n", model_suptype)
	if (strlen(inflation_line) > 0) {
		printf("%s\n", inflation_line)
	}
	printf("Regime switching:        %s  \n", switching_type)
	printf("Number of observations = %9.0f \n", model.n)
	printf("Log likelihood         = %9.4f \n", model.logLik)
	printf("McFadden pseudo R2     = %9.4f \n", model.R2)
	printf("LR chi2(%2.0f)            = %9.4f \n", model.df - model.df_null, model.chi2)
	printf("Prob > chi2            = %9.4f \n", model.chi2_pvalue)
	printf("AIC                    = %9.4f \n" , model.AIC)
	printf("BIC                    = %9.4f \n" , model.BIC)
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
	
	model = postDescribeModel(model, robust, allvars, corresp)
	
	model.eqnames = J(1, cols(tokens(xnames)) + 2, "Regime equation"),  J(1, cols(tokens(zpnames)) + model.ncatp, "Outcome equation (+)"),  J(1, cols(tokens(znnames)) + model.ncatn, "Outcome equation (-)")
	model.parnames = tokens(xnames), "/cut1", "/cut2", tokens(zpnames),  "/cut" :+ strofreal(1..model.ncatp), tokens(znnames),  "/cut" :+ strofreal(1..model.ncatn)
	
	if (correlated) {
		model.eqnames = model.eqnames, J(1, 2, "Correlation coefficients")
		model.parnames = model.parnames, "rho(+)", "rho(-)"
	}
	
	passModelToStata(model)
	
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
	
	model = postDescribeModel(model, robust, allvars, corresp)
	
	model.eqnames = J(1, cols(tokens(xnames))+2, "Regime equation"),  J(1, cols(tokens(zpnames)) + model.ncatn - 1, "Outcome equation (+)"),  J(1, cols(tokens(znnames)) + model.ncatn - 1, "Outcome equation (-)")
	model.parnames = tokens(xnames), "/cut1", "/cut2", tokens(zpnames),  "/cut" :+ strofreal(1..(model.ncatp-1)), tokens(znnames),  "/cut" :+ strofreal(1..(model.ncatn-1))
	
	if (correlated) {
		model.eqnames = model.eqnames, J(1, 2, "Correlation coefficients")
		model.parnames = model.parnames, "rho(+)", "rho(-)"
	}
	
	passModelToStata(model)
	
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
	model.yname = yname
	model.xnames = xnames
	model.znames = znames
	
	model = postDescribeModel(model, robust, allvars, corresp)
	
	model.eqnames = J(1, cols(tokens(xnames)) + 1, "Regime equation"), J(1, cols(tokens(znames)) + rows(model.allcat)-1, "Outcome equation")
	model.parnames = tokens(xnames), "/cut1", tokens(znames), "/cut" :+ strofreal(1..(rows(model.allcat)-1))
	
	if (correlated) {
		model.eqnames = model.eqnames, "Correlation coefficient"
		model.parnames = model.parnames, "rho"
	}
		
	passModelToStata(model)
	
	return(model)
}

function escape_stripes(stripes) {
	// workaround: stata does not allow colstripes containing dots
	colstripes = subinstr(stripes, "=.", "=0.")
	colstripes = subinstr(colstripes, "=-.", "=-0.")
	colstripes = subinstr(colstripes, ".", ",")
	return(colstripes)
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
	return(colstripes)
}

function output_matrix(matrix_name, matrix_value, rowstripes, colstripes){
	rowstripes_new = escape_stripes(rowstripes)
	colstripes_new = escape_stripes(colstripes)
	st_matrix(matrix_name, matrix_value)
	st_matrixrowstripe(matrix_name, (J(rows(rowstripes_new), 1, ""), rowstripes_new))
	st_matrixcolstripe(matrix_name, (J(rows(colstripes_new), 1, ""), colstripes_new))
}

function output_mesetp(me, se, rowstripes, colstripes) {
	t = me :/ se
	pval = (1:-normal(abs(t))) :* 2
	output_matrix("me",     me, rowstripes, colstripes)
	output_matrix("se",     se, rowstripes, colstripes)
	output_matrix("t",       t, rowstripes, colstripes)
	output_matrix("pval", pval, rowstripes, colstripes)
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

function CNOPmargins(class CNOPModel scalar model, string atVarlist, zeroes, regime) {
	xzbar = model.XZmedians
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
	output_matrix("at", xzbar, " ", model.XZnames')
	
	rowstripes = model.XZnames'
	colstripes = get_colstripes(model.model_class, loop, model.allcat, model.infcat)
	
	mese = generalMEwithSE(xzbar, model, loop)
	kxz = cols(xzbar)
	me = mese[1::kxz,]
	se = mese[(1::kxz) :+ kxz,]
	
	output_mesetp(me, se, rowstripes, colstripes)
	
	// now the printing part! 
	"Evaluated at:"
	print_matrix(xzbar, ., model.XZnames)
	""
	if (zeroes) {
		"Marginal effects of all variables on the probabilities of different types of zeros"
	} 
	else if (regime) {
		"Marginal effects of all variables on the probabilities of different latent regimes"
	}
	else {
		"Marginal effects of all variables on the probabilities of different outcomes"
	}
	print_matrix(me, rowstripes, colstripes)
	""
	"Standard errors of marginal effects"
	print_matrix(se, rowstripes, colstripes)
}


function CNOPprobabilities(class CNOPModel scalar model, string atVarlist, zeroes, regime) {
	xz_from = model.XZmedians
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
	
	output_matrix("at", xz_from, " ", model.XZnames')

	colstripes = get_colstripes(model.model_class, loop, model.allcat, model.infcat)
	rowstripes = " " // rowstripes made invisible
	mese = generalPredictWithSE(xz_from, model, loop)
	me = mese[1,]
	se = mese[2,]
	output_mesetp(me, se, rowstripes, colstripes)
	
	// now the printing part! 
	"Evaluated at:"
	print_matrix(xz_from, ., model.XZnames)
	""
	if (zeroes) {
		"Predicted probabilities of different types of zeros"
	} 
	else if (regime) {
		"Predicted probabilities of different latent regimes"
	}
	else {
		"Predicted probabilities of different outcomes"
	}
	print_matrix(me, ., colstripes)
	""
	"Standard errors of the probabilities"
	print_matrix(se, ., colstripes)
}

function CNOPcontrasts(class CNOPModel scalar model, string atVarlist, string toVarlist, zeroes, regime) {
	xz_from = model.XZmedians
	xz_to = model.XZmedians
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
	
	output_matrix("between", xz_from \ xz_to, "from" \ "to", model.XZnames')
	
	colstripes = get_colstripes(model.model_class, loop, model.allcat, model.infcat)
	rowstripes = " " // rowstripes made invisible
	
	mese = generalContrastsWithSE(xz_from, xz_to, model, loop)
	me = mese[1,]
	se = mese[2,]
	output_mesetp(me, se, rowstripes, colstripes)
	
	// now the printing part! 
	"Evaluated between"
	print_matrix(xz_from \ xz_to, "from" \ "to", model.XZnames')
	""
	if (zeroes) {
		"Contrasts of the predicted probabilities of different types of zeros"
	} 
	else if (regime) {
		"Contrasts of the predicted probabilities of different latent regimes"
	}
	else {
		"Contrasts of the predicted probabilities of different outcomes"
	}
	print_matrix(me, ., colstripes)
	""
	"Standard errors of the contrasts"
	print_matrix(se, ., colstripes)
	
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
	if (output == "mode" | output == "choice") {
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

void vuong_calc(| ll_diff, k_1, k_2){
	if (rows(ll_diff) < 1 || ll_diff == .) {
		ll_diff = st_matrix("ll_diff")
		k_1 = strtoreal(st_local("k_1"))
		k_2 = strtoreal(st_local("k_2"))
	}
	mean_diff = mean(ll_diff)
	std_diff = sqrt(variance(ll_diff))
	n_obs = rows(ll_diff)
	vuong = mean_diff / (std_diff / sqrt(n_obs))
	
	// AIC and BIC corrections
	vuongAIC = (mean_diff - (k_1-k_2) / n_obs) / (std_diff / sqrt(n_obs))
	vuongBIC = (mean_diff - (k_1-k_2) * log(n_obs) / (2 * n_obs)) / (std_diff / sqrt(n_obs))
	
	pvalue = 1-normal(vuong)
	pvalueAIC = 1-normal(vuongAIC)
	pvalueBIC = 1-normal(vuongBIC)
	
	sprintf("Mean difference in log likelihood                  %9.4f", mean_diff)
	sprintf("Standard deviation of difference in log likelihood %9.4f", std_diff)
	sprintf("Number of observations                             %9.0f", n_obs)
	sprintf("Vuong test statistic                           z = %9.4f", vuong)
	sprintf("P-Value                                     Pr>z = %9.4f", pvalue)
	sprintf("   with AIC (Akaike) correction                z = %9.4f", vuongAIC)
	sprintf("P-Value                                     Pr>z = %9.4f", pvalueAIC)
	sprintf("   with BIC (Schwarz) correction               z = %9.4f", vuongBIC)
	sprintf("P-Value                                     Pr>z = %9.4f", pvalueBIC)
	
	st_numscalar("mean_diff", mean_diff)
	st_numscalar("std_diff", std_diff)
	st_numscalar("n_obs", n_obs)
	st_numscalar("vuong", vuong)
	st_numscalar("vuongAIC", vuongAIC)
	st_numscalar("vuongBIC", vuongBIC)
	st_numscalar("pvalue", pvalue)
	st_numscalar("pvalueAIC", pvalueAIC)
	st_numscalar("pvalueBIC", pvalueBIC)
}


void classification_calc_large(fact_varname, pred_varname, touse) {
	// access Stata variables
	st_view(fact = ., ., fact_varname, touse)
	st_view(pred = ., ., pred_varname, touse)
	
	// construct the classification table
	all_levels = uniqrows(fact \ pred)
	ncat = rows(all_levels)
	confmat = J(ncat, ncat, 0)
	for (i=1; i <= ncat; i++) {
		for (j=1; j <= ncat; j++) {
			confmat[i, j] = sum((fact :== all_levels[i]) :& (pred :== all_levels[j]))
		}
	}
	// todo: calculate margins of the table, add row and column titles
	//print_matrix(confmat, strofreal(all_levels), strofreal(all_levels), ., ., ., 0)
	st_matrix("cells", confmat)
	st_matrix("labels", all_levels)
	classification_calc("cells", "labels", "result")
}

void classification_calc(cells_matname, labels_matname, result_matname) {
	// analyze classification table
	
	cells = st_matrix(cells_matname)
	labels = st_matrix(labels_matname)
	
	total = sum(cells)
	
	ap = rowsum(cells)
	an = total :- ap
	pp = colsum(cells)'
	pn = total :- pp
	
	tp = diagonal(cells)
	fp = pp - tp
	fn = ap - tp
	tn = pn - fn
	
	noise = fp :/ an
	recall = tp :/ ap
	precision = tp :/ pp
	n2s = noise :/ recall
	
	result = precision, recall, n2s
	colnames = "Precision" \  "Recall" \  "Adjusted noise-to-signal ratio"
	rownames = strofreal(labels)
	
	output_matrix(result_matname, result, rownames, colnames) 
	
	rowtitle = ("Actual" \ "outcomes")
	print_matrix(result, strofreal(labels), colnames, uline=., lline=., mline=., digits=., rowtitle=rowtitle)
}

void print_matrix(contents, rownames, colnames, | uline, lline, mline, digits, rowtitle, coltitle) {
	// because Stata cannot display matrices with dots in colnames, we need our own printing function!
	n = rows(contents)
	m = cols(contents)
	if (rownames == . | rows(rownames) == 0) {
		rowname_width = 0
		rowname_flag = 0
	} else {
		rowname_width = max(strlen(rownames) \ 10)
		rowname_flag = 1
	}
	if (uline == . | rows(uline) == 0) {
		uline = 0
	} 
	if (lline == . | rows(lline) == 0) {
		lline = 0
	}
	if (mline == . | rows(mline) == 0) {
		mline = (n > 1)
	}
	if (digits == . | rows(digits) == 0) {
		digits = 4
	}
	_colnames = colnames
	if (cols(_colnames) > 1){
		_colnames = _colnames'
	}
	
	if (rowtitle == . | rows(rowtitle) == 0) {
		rowtitle_rows = 0
	} else {
		// todo: ensure that rowname_flag is true
		rowtitle_rows = rows(rowtitle)
		rowname_width = max((strlen(rowtitle) \ rowname_width))
	}
	if (coltitle == . | rows(coltitle) == 0) {
		coltitle_rows = 0
	} else {
		// todo: ensure that rowname_flag is true
		coltitle_rows = rows(coltitle)
	}
	
	colwidths = rowmax((strlen(_colnames) :+ 3 , J(rows(_colnames), 1, 6)))
	// todo: support word wrap for long colnames and maybe row and col titles
	// todo: make colwidths depend on the contents
	// todo: support lines before totals
	numberf = strofreal(digits) + "f"
	if (rowname_flag) {
		hline = "{hline " + strofreal(rowname_width+1)+ "}{c +}{hline " + strofreal(sum(colwidths :+ 1) + 2)+ "}\n"
	} else {
		hline = "{hline " + strofreal(rowname_width+1+1+sum(colwidths :+ 1) + 2) + "}\n"
	}
	// print header
	if (uline) {
		printf(hline)
	}
	
	if (rowtitle_rows > 1) {
		for(i=1; i <= rowtitle_rows; i++) {
			// todo: take into accoutn possible difference in vlines
			printf("%" + strofreal(rowname_width) + "s {c |}", rowtitle[i])
			// todo: make coltitle centered
			if (coltitle_rows > 0) {
				coltitle_current =  i + coltitle_rows - rowtitle_rows + 1
				if ((coltitle_current > 0) & (coltitle_current <= coltitle_rows)) {
					printf(coltitle[coltitle_current])
				}
			}
			if (i < rowtitle_rows) {
				printf("\n")
			}
		}
	} else if (rowname_flag) {
		printf("%" + strofreal(rowname_width) + "s {c |} ", "")
	}
	for(j=1; j<=m; j++){
		printf("%" + strofreal(colwidths[j]) + "s ", colnames[j])
	}
	printf("\n")
	if (mline) {
		printf(hline)
	}
	// print the rest of the table
	for(i=1; i<=n; i++) {
		if (rowname_flag) {
			printf("%" + strofreal(rowname_width)+ "s {c |} ", rownames[i])
		}
		for(j=1; j<=m; j++){
			printf("%" + strofreal(colwidths[j]) + "." + numberf + " ", contents[i, j])
		}
		printf("\n")
	}
	if (lline) {
		printf(hline)
	}
}

end

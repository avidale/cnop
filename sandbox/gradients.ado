version 12

mata

/*
a bit of theory:

if f(a,b) = Normal(a, b | r), then

df/da = normalden(a) * cumnormal((b-ar)/(1-r^2))

*/

/*
INCOMPLETE: insert these constants into each optimization
*/
CNOP_GLOBAL_CONST_MAXITER = 30
CNOP_GLOBAL_CONST_PTOL = 1e-6
CNOP_GLOBAL_CONST_VTOL = 1e-7
CNOP_GLOBAL_CONST_NRTOL = 1e-5
CNOP_GLOBAL_CONST_LAMBDA = 0.0



////////////////////////////////////////////////////////////////
// support functions
////////////////////////////////////////////////////////////////

function punishSort(x, | delta, diff) {
	// the function returns vector x, if it is sorted in ascending order
	// otherwise, it returns modified x in incrieasing order with very small distances between swapped points
	// the goal is to avoid negative probability
	if(args()==1){
		delta = 10^-5
	}
	mu = x;
	n = length(mu);
	diff  = J(n, 1, 0)
	i = 1;
	while(i<n){
		if(mu[i+1] <= mu[i]) {
			diff[i+1] = mu[i+1] - mu[i] + delta
			mu[i+1] = mu[i] + delta
		}
		i=i+1;
	}
	s = sum(diff:^2)
	if (s>0) {
		//x'
		//s
	}
	return(mu);
}

function codeIncreasingSequence(decoded_sequence) {
	/* incomplete: check if increasing! */
	coded_sequence = decoded_sequence
	n = length(coded_sequence)
	for (i = 2; i <= n; i++) {
		coded_sequence[i] = log(decoded_sequence[i] - decoded_sequence[i-1])
	}
	return(coded_sequence)
}
function decodeIncreasingSequence(coded_sequence) {
	decoded_sequence = coded_sequence
	n = length(coded_sequence)
	for (i = 2; i <= n; i++) {
		next = coded_sequence[i]
		/* for numerical stability we have to bound exponentials */
		next = min((700, max((-700, next))))
		decoded_sequence[i] = decoded_sequence[i-1] + exp(next)
	}
	return(decoded_sequence)
}
function codeIncreasingSequenceJacobian(coded_sequence) {
	n = length(coded_sequence)
	expx = exp(coded_sequence)
	jacobian = J(n, n, 0)
	for (j = 1; j <= n; j++) {
		jacobian[1, j] = 1;
		for (i = 2; i <= j; i++) {
			jacobian[i, j] = expx[i]
		}
	}
	return(jacobian)
}

// get argsort of the nonzero positions
function selectOrderedIndex(positions) {
	// find unordered nonzero positions
	where = selectindex(positions)
	index = where[order(positions[where]', 1)']
	return(index)
}


function coeffOP(x, ycateg, ncat,| quiet, startbmu, lambda, maxiter, ptol, vtol, nrtol) {
	ntry = 0
	if (args() < 4 || rows(quiet) < 1 || quiet == .) {
		quiet = 0
	}
	if (args() <5 || rows(startbmu) < 1 || startbmu == .) {
		ntry = 1
		cumprob 	= runningsum(mean(ycateg));
		start_mu1	= invnormal(cumprob[1::ncat-1])'; // assume intercept as in a model with no slope
		start_b1	= invsym(x'*x)*x'*(ycateg * invnormal(0.5 :* cumprob + 0.5 :* (0, cumprob[1::ncat-1]))'); // try to predict median for each category 
		startbmu 	= start_b1 \ start_mu1;
	}
	if (args() < 6 || rows(lambda) < 1 || lambda ==.) {
		external CNOP_GLOBAL_CONST_LAMBDA
		lambda = CNOP_GLOBAL_CONST_LAMBDA
	}
	if (args() < 7 || rows(maxiter) < 1 || maxiter ==.) {
		external CNOP_GLOBAL_CONST_MAXITER
		maxiter = CNOP_GLOBAL_CONST_MAXITER
	}
	if (args() < 8 || rows(ptol) < 1 || ptol ==.) {
		external CNOP_GLOBAL_CONST_PTOL
		ptol = CNOP_GLOBAL_CONST_PTOL
	}
	if (args() < 9 || rows(vtol) < 1 || vtol ==.) {
		external CNOP_GLOBAL_CONST_VTOL
		vtol = CNOP_GLOBAL_CONST_VTOL
	}
	if (args() < 10 || rows(nrtol) < 1 || nrtol ==.) {
		external CNOP_GLOBAL_CONST_NRTOL
		nrtol = CNOP_GLOBAL_CONST_NRTOL
	}

	S = optimize_init()
	if(quiet){
		optimize_init_tracelevel(S, "none")
		optimize_init_verbose(S, 0)
	}
	optimize_init_argument(S, 1, x)
	optimize_init_argument(S, 2, ycateg)
	optimize_init_argument(S, 3, ncat)
	optimize_init_argument(S, 4, lambda)
	optimize_init_evaluator(S, &_op_optim())
	optimize_init_evaluatortype(S, "gf0") // unresolved errors=(
	optimize_init_params(S, (startbmu'))
	
	optimize_init_conv_maxiter(S, maxiter)
	optimize_init_conv_ptol(S, ptol)
	optimize_init_conv_vtol(S, vtol)
	optimize_init_conv_nrtol(S, nrtol)
	
	optimize_init_conv_warning(S, "off")
	
	//optimize_init_singularHmethod(S, "hybrid")
	optimize_init_verbose(S, 0)
	code = _optimize(S)
	
	while (code != 0) {
		"OP: smth happened"
		optimize_result_errortext(S)
		cumprob 	= runningsum(mean(ycateg));
		ycateg
		cumprob
		start_mu1	= invnormal(cumprob[1::ncat-1])'; // assume intercept as in a model with no slope
		start_mu1
		start_b1	= invsym(x'*x)*x'*(ycateg * invnormal(0.5 :* cumprob + 0.5 :* (0, cumprob[1::ncat-1]))'); // try to predict median for each category 
		startbmu 	= start_b1 \ start_mu1;
		start_b
		if(ntry == 0) {
			// use the best seeming initial values
		} else if(ntry == 1) {
			start_b1 = start_b1 * 0
			startbmu 	= start_b1 \ start_mu1
			optimize_init_params(S, (startbmu'))
		} else if(ntry==2) {
			optimize_init_evaluatortype(S, "gf0")
		} else if(ntry==3) {
			optimize_init_technique("bfgs")
		} else {
			"CNOP ordered probit estimation: no way to converge"
			break;
		}
		ntry = ntry+1
		code = _optimize(S)'
	}
	
	answer = optimize_result_params(S)'
	
	return(answer)
}

function coeffOP2(x, ycateg, ncat,| startbmu) {
	ntry = 0

	if (args() == 3) {
		ntry = 1
		cumprob = runningsum(mean(ycateg));
		mu		= invnormal(cumprob[1::ncat-1])';
		nu 		= mu[1] \ log(mu[2::ncat-1] - mu[1::ncat-2])
		b		= invsym(x'*x)*x'*(ycateg * invnormal(0.5 :* cumprob + 0.5 :* (0, cumprob[1::ncat-1]))'); 
		startbnu 	= b \ nu;
	}

	S = optimize_init()
	optimize_init_argument(S, 1, x)
	optimize_init_argument(S, 2, ycateg)
	optimize_init_argument(S, 3, ncat) 
	optimize_init_evaluator(S, &_op_optim_expnu())
	optimize_init_evaluatortype(S, "gf1debug")
	optimize_init_params(S, (startbnu'))
	
	optimize_init_conv_ptol(S, 10^-8)
	optimize_init_conv_vtol(S, 10^-9)
	optimize_init_conv_nrtol(S, 10^-7)
	
	optimize_init_conv_warning(S, "off")
	
	//optimize_init_singularHmethod(S, "hybrid")
	optimize_init_verbose(S,0)
	code = _optimize(S)
	//code
	
	while (code != 0) {
		"smth happened"
		optimize_result_errortext(S)
		cumprob 	= runningsum(mean(ycateg));
		ycateg
		cumprob
		start_mu1	= invnormal(cumprob[1::ncat-1])'; // assume intercept as in a model with no slope
		start_mu1
		start_b1	= invsym(x'*x)*x'*(ycateg * invnormal(0.5 :* cumprob + 0.5 :* (0, cumprob[1::ncat-1]))'); // try to predict median for each category 
		startbmu 	= start_b1 \ start_mu1;
		start_b
		if(ntry == 0) {
			// use the best seeming initial values
		} else if(ntry == 1) {
			start_b1 = start_b1 * 0
			startbmu 	= start_b1 \ start_mu1
			optimize_init_params(S, (startbmu'))
		} else if(ntry==2) {
			optimize_init_evaluatortype(S, "gf0")
		} else if(ntry==3) {
			optimize_init_technique("bfgs")
		} else {
			"CNOP ordered probit estimation: no way to converge"
			break;
		}
		ntry = ntry+1
		code = _optimize(S)'
	}
	
	answer = optimize_result_params(S)'
	
	return(answer)
}

void _op_optim_expnu(todo, params, x, q, ncat, v, g, H) {
	// nu = mu[1] \ log(mu[2::ncat-1] - mu[1::ncat-2])

	k 		= cols(x)
	n		= rows(x)
	gama 	= params[1::k]'
	nu 		= params[(k+1)::length(params)]
	dif 	= exp(nu)'
	dif[1]	= nu[1]
	mu		= runningsum(dif)
	xg 		= x * gama
	prob 	= normal(mu'[J(n,1,1),] :- xg) , J(n,1,1)
	prob[ , 2::ncat] = prob[ , 2::ncat] - prob[ , 1::(ncat-1)]
	v		= rowsum(log(prob) :* q) 
	if (todo == 1) {
		g = J(n, length(params), 0)
		dprob = normalden(mu'[J(n,1,1),] :- xg)
		g[, 1::k] = rowsum(((dprob, J(n, 1, 0)) - (J(n, 1, 0), dprob)) :* q) :* x
		dif[1] = 1
		g[, k+1::k+ncat-1] = (dprob :* (q[,1::ncat-1] - q[,2::ncat])) * (lowertriangle(J(ncat-1, ncat-1, 1)) :* dif')
		g = g :/ rowsum(prob :* q) 
	}
}

function tokenCoincidence(string1, string2){
	/*
	example:
	tokenCoincidence("a b c", "d e f")
	1 2 3 0 0 0
	0 0 0 1 2 3
	tokenCoincidence("a b c", "d b a")
	1 2 3 0
	3 2 0 1
	can be used with st_varindex and tokens
	st_matrix, st_matrixcolstripe for results
	*/
}
function calculateGradient(func, x, otherargs) {
	dx = min(x * 1e-08, 1e-10)
	return( ((*func)(x+dx)-(*func)(x)) / dx )
}


function referenceCategory(currentCategory, uniqueCategories) {
	// the function intended for ordered categorical variables. For ordinary dummies works by default.
	ncat = length(uniqueCategories);
	catIndex = selectindex(uniqueCategories == currentCategory);
	if (catIndex == ncat) {
		return (uniqueCategories[catIndex - 1]);
	} else {
		return (uniqueCategories[catIndex + 1]);
	}
}
function referenceValues(xzbar, xzSet, dummies) {
	refValues = xzbar;
	for (i = 1; i <= cols(xzbar); i++) {
		if (dummies[i] != 0) {
			if (dummies[i] == 1) {
				uniqueCategories = uniqrows(xzSet[,i]);
			} else if (dummies[i] == 2) {
				uniqueCategories = range(min(xzSet[,i]), max(xzSet[,i]), 0.25);
			} 
			refValues[i] = referenceCategory(xzbar[i], uniqueCategories);
		}
	}
}



////////////////////////////////////////////////////////////////
// general functions
////////////////////////////////////////////////////////////////
function generalPredict(real matrix xzvalues, class CNOPModel scalar model, real scalar loop) {
	// loop: 1 - probabilities of y labels, 2 - decompositions of zeros, 3 - probabilities of -1 - 0 -1 regimes 
	if (sum(model.model_class :== ("CNOP", "CNOPC", "NOP", "NOPC")) > 0){
		x	= xzvalues[,selectOrderedIndex(model.corresp[1,])]
		zp	= xzvalues[,selectOrderedIndex(model.corresp[2,])]
		zn	= xzvalues[,selectOrderedIndex(model.corresp[3,])]
	} else if(sum(model.model_class :== ("MIOPR", "MIOPRC")) > 0){
		x	= xzvalues[,selectOrderedIndex(model.corresp[1,])]
		z	= xzvalues[,selectOrderedIndex(model.corresp[2,])]
	}
	if (model.model_class == "CNOP") {
		prediction = MLcnop(model.params, x, zp, zn, ., model.ncat, model.infcat, loop)
	} else if(model.model_class == "CNOPC") {
		prediction = MLcnopc(model.params, x, zp, zn, ., model.ncat, model.infcat, loop)
	} else if(model.model_class == "MIOPR"){
		prediction = MLmiopr(model.params, x, z, ., model.ncat, model.infcat, loop)
	} else if(model.model_class == "MIOPRC"){
		prediction = MLmioprc(model.params, x, z, ., model.ncat, model.infcat, loop)
	} else if (model.model_class == "NOP") {
		prediction = MLnop(model.params, x, zp, zn, ., model.ncat, model.infcat, loop)
	} else if (model.model_class == "NOPC") {
		prediction = MLnopc(model.params, x, zp, zn, ., model.ncat, model.infcat, loop)
	} else if (model.model_class == "OP") {
		prediction = MLop(model.params, xzvalues, ., model.ncat, loop)
	} else {
		"predict() not determined for model " + model.model_class
		prediction = .
	}
	return(prediction)
}

function generalPredictWrapper(real matrix params, real matrix xzvalues, class CNOPModel scalar model, real scalar loop, real matrix returnedPred)	{
	real_params = model.params * 1
	model.params = params'
	
	returnedPred = generalPredict(xzvalues, model, loop)
	model.params = (real_params * 1)
	
	returnedPred = rowshape(returnedPred, 1)
}

function generalPredictWithSE(real matrix xzvalues, class CNOPModel scalar model, real scalar loop) {
	generalPredictWrapper(model.params', xzvalues, model, loop, probs = .)
	probs = rowshape(probs, 1)
	
	nc = cols(probs)
	nr = rows(probs)
	D = deriv_init()
	deriv_init_evaluator(D, &generalPredictWrapper())
	deriv_init_evaluatortype(D, "t")
	
	deriv_init_params(D, model.params')
	deriv_init_argument(D, 1, xzvalues)
	deriv_init_argument(D, 2, model)
	deriv_init_argument(D, 3, loop)
	
	errorcode = _deriv(D, 1)
	if (errorcode != 0) {
		" Error in numerical PREDICT (SE) differentiation " + strofreal(errorcode)
		// INCOMPLETE: do something if differentiation fails; research whether it can be corrected
	}
	/*grad = rowshape(deriv_result_Jacobian(D), np)*/
	grad = deriv_result_Jacobian(D)'
	
	se = J(nr,nc,.)
	
	if (model.robust == 1) {
		V = model.V_rob
	} else {
		V = model.V
	}
	
	for(i=1; i<=nr;++i) {
		for(j=1; j<=nc;++j){
			gradrow= grad[,(i-1)*nc + j]'
			se[i,j] = gradrow * V * gradrow'
		}
	}
	se = sqrt(se)
	
	return (probs \ se)
		
}

//
function generalME(real matrix params, real vector xzbar, class CNOPModel scalar model, real vector dummies, real scalar loop, real matrix _returnedME) {
	// xzbar, pmf, params, ncat, infcat, dummies, corresp, rawME, _returnedME
	// wanna have derivative of PMF!!!

	if (model.model_class == "CNOP") {
		ME = cnop_me_raw(params', xzbar, model.ncat, model.infcat, model.corresp, loop)
	} else if(model.model_class == "CNOPC") {
		ME = cnopc_me_raw(params', xzbar, model.ncat, model.infcat, model.corresp, loop)
	} else if(model.model_class == "MIOPR") {
		ME = miopr_me_raw(params', xzbar, model.ncat, model.infcat, model.corresp, loop)
	} else if(model.model_class == "MIOPRC") {
		ME = mioprc_me_raw(params', xzbar, model.ncat, model.infcat, model.corresp, loop)
	} else if(model.model_class == "NOP") {
		ME = nop_me_raw(params', xzbar, model.ncat, model.infcat, model.corresp, loop)
	} else if(model.model_class == "NOPC") {
		ME = nopc_me_raw(params', xzbar, model.ncat, model.infcat, model.corresp, loop)
	} else if(model.model_class == "OP") {
		ME = op_me_raw(params', xzbar, model.ncat, model.infcat, model.corresp, loop)
	} else {
		"ME not determined for model of class: " + model.model_class
		ME = .
	}
	params_tmp = model.params
	model.params = params'
	
	probBase = generalPredict(xzbar, model, loop)
	for (i = 1; i <= cols(xzbar); i++) {
		if (dummies[i] != 0) {
			xzbarOther = xzbar
			xzbarOther[i] = 1 - xzbarOther[i]
			probOther = generalPredict(xzbarOther, model, loop)
			ME[i, ] =  probOther - probBase
		}
	}
	model.params = params_tmp
	_returnedME = rowshape(ME, 1); // change the by-reference argument!!!
}
//
function generalMEwithSE(real vector xzbar, class CNOPModel scalar model, real vector dummies, real scalar loop) {
	
	generalME(model.params', xzbar, model, dummies, loop, me=.) // the function changes the last argument
	me = rowshape(me, length(xzbar))
	nc = cols(me)
	nr = rows(me)
	np = rows(model.params)
	D = deriv_init()
	deriv_init_evaluator(D, &generalME())
	deriv_init_evaluatortype(D, "t")
	
	deriv_init_params(D, model.params')
	deriv_init_argument(D, 1, xzbar)
	deriv_init_argument(D, 2, model)
	deriv_init_argument(D, 3, dummies)
	deriv_init_argument(D, 4, loop)
	
	errorcode = _deriv(D, 1)
	if (errorcode != 0) {
		" Error in numerical ME differentiation " + strofreal(errorcode)
		model.params
		// INCOMPLETE: do something if differentiation fails; research whether it can be corrected
	}
	/*grad = rowshape(deriv_result_Jacobian(D), np)*/
	grad = deriv_result_Jacobian(D)'
	
	se = J(nr,nc,.)
	
	if (model.robust == 1) {
		V = model.V_rob
	} else {
		V = model.V
	}
	for(i=1; i<=nr;++i) {
		for(j=1; j<=nc;++j){
			gradrow= grad[,(i-1)*nc + j]'
			se[i,j] = gradrow * V * gradrow'
		}
	}
	se = sqrt(se)
	model.me = me
	model.mese = se
	
	return (me \ se)
}



function generalContrasts(real matrix params, real vector xz_from, real vector xz_to, class CNOPModel scalar model, real scalar loop, real matrix _returnedME) {
	// I make tricks with parameters because I want the prediction to depend on the parameters vector to make delta-metod possible
	// returned value is the row of differences between y1...yn probabilities at two xz points
	real_params = model.params'
	model.params = params'
	probBase = generalPredict(xz_from, model, loop)
	probOther = generalPredict(xz_to, model, loop)
	model.params = real_params'
	_returnedME = rowshape(probOther - probBase, 1); // change the by-reference argument!!!
}

function generalContrastsWithSE(real vector xz_from, real vector xz_to, class CNOPModel scalar model, real scalar loop) {
	// here ME is just a row - difference between two probability vectors
	generalContrasts(model.params',  xz_from, xz_to, model, loop, me=.)
	
	D = deriv_init()
	deriv_init_evaluator(D, &generalContrasts())
	deriv_init_evaluatortype(D, "t")
	
	deriv_init_params(D, model.params')
	deriv_init_argument(D, 1, xz_from)
	deriv_init_argument(D, 2, xz_to)
	deriv_init_argument(D, 3, model)
	deriv_init_argument(D, 4, loop)
	
	me = rowshape(me, 1)
	nc = cols(me)
	nr = rows(me)
	np = rows(model.params)
	grad = rowshape(deriv(D, 1), np)
	
	se = J(nr,nc,.)
	
	if (model.robust == 1) {
		V = model.V_rob
	} else {
		V = model.V
	}
	
	for(i=1; i<=nr;++i) {
		for(j=1; j<=nc;++j){
			gradrow= grad[,(i-1)*nc + j]'
			se[i,j] = gradrow * V * gradrow'
		}
	}
	se = sqrt(se)
	
	return (me \ se)
	
}
////////////////////////////////////////////////////////////////
// ORDERED PROBIT MODEL
////////////////////////////////////////////////////////////////
// likelihood of a OP model
function MLop(params, x, q, ncat, | loop) {
	k 		= cols(x)
	n		= rows(x)
	gama 	= params[1::k]
	mu 		= params[(k+1)::length(params)] 	
	mu		= punishSort(mu)
	xg 		= x * gama
	/*
	prob	= J(rows(x), ncat, 0)
	prob[.,1] 	= normal(mu[1] :- xg)
	for(i = 2; i<=ncat-1;i++){
		prob[.,i]	 = normal(mu[i] :- xg) - normal(mu[i-1] :- xg)
	}
	prob[.,ncat] = 1:-normal(mu[ncat-1] :- xg)
	*/
	// more modern: create F(mu-xg) as a matrix and subtract from itself
	prob = normal(mu'[J(n,1,1),] :- xg) , J(n,1,1)
	prob[ , 2::ncat] = prob[ , 2::ncat] - prob[ , 1::(ncat-1)]
	
	if(loop==1){
		return(prob)
	}else{
		col_logl	= rowsum(log(prob) :* q) //added!
		logl 	= colsum(col_logl)
		return(col_logl) // added!
	}
}

// gradient of likelihood of a OP model
function Jacop(params, x, q, ncat) {
	k 		= cols(x)
	n 		= rows(x)
	gama 	= params[1::k]
	mu 		= params[(k+1)::length(params)] 	// must be ordered!!!!!!!!!
	xg 		= x * gama
	
	// the most horrible loops I have ever seen
	grad 	= (-normalden(mu[1]:-xg) :* x :* q[.,1] :/ normal( mu[1]:-xg ) ) , ( normalden(mu[1]:-xg) :* q[.,1] :/ normal( mu[1]:-xg ) )
	
	if(ncat>2){
		grad = (grad), (J(n, ncat-2,0))
	}
	for( j=2; j <= ncat - 1; j++){
		g_ij = -(normalden(mu[j] :- xg) :- normalden(mu[j-1] :- xg)) :* x :* q[.,j] :/(normal(mu[j] :- xg) :- normal(mu[j-1] :- xg));
		if(j == 2){
			g_ij = g_ij, (-normalden(mu[j-1]:-xg) :* q[.,j] :/ (normal(mu[j] :- xg) - normal(mu[j-1] :- xg))),  (normalden(mu[j]:-xg) :* q[.,j] :/ (normal(mu[j] :- xg) - normal(mu[j-1] :- xg)))
			if( ncat > 3){
				 g_ij = g_ij , J(n,ncat-1-2,0)
			}
		} else if(j < ncat - 1){
			g_ij	= g_ij, J(n, j-2, 0),  (-normalden(mu[j-1]:-xg) :* q[.,j] :/ (normal(mu[j] :- xg) - normal(mu[j-1] :- xg))), (normalden(mu[j]:-xg) :* q[.,j] :/ (normal(mu[j] :- xg) - normal(mu[j-1] :- xg))), J(n, ncat-1-j, 0)
		} else {	
			g_ij	= g_ij, J(n, j-2, 0),  (-normalden(mu[j-1]:-xg) :* q[.,j] :/ (normal(mu[j] :- xg) - normal(mu[j-1] :- xg))), (normalden(mu[j]:-xg) :* q[.,j] :/ (normal(mu[j] :- xg) - normal(mu[j-1] :- xg)))
		}
		grad = grad + g_ij
	}
	
	g_ij = normalden(mu[ncat-1] :- xg) :* q[., ncat] :* x :/ (1 :- normal(mu[ncat-1] :- xg));
	if(ncat > 2){
		g_ij = g_ij, J(n, ncat - 2, 0)
	}
	g_ij = g_ij , ( -normalden(mu[ncat-1] :- xg) :* q[.,ncat] :/ (1 :- normal(mu[ncat-1] :- xg))) 
	
	grad = grad + g_ij;
	
	return(grad);
}

// ME for an OP model
void _op_optim(todo, params, x, q, ncat, lambda, v, g, H) {
	v		= MLop(params', x, q, ncat, 0)
	n = rows(q)
	v = v :- lambda / n * sum(params:^2) 
	if(todo==1){
		g 	= Jacop(params', x, q, ncat)
		g = g :- 2 * lambda * params /* INCOMPLETE: do I need to transpose params? */ 
	}
}
function op_me(params, xbar, q, ncat, | loop) {	 // loop as a parameter
	k 		= cols(xbar)
	gama 	= params[1::k]
	mu 		= params[(k+1)::length(params)] 	// must be ordered!!!!!!!!!
	xg 		= xbar * gama
	
	mes		= normalden(mu[1] :- xg) :* -gama;	// ME for y=y0 with respect to all x' : k \times 1 vector
	
	for (j=2; j<=ncat-1; j++){
		mes = (mes) , ((normalden(mu[j] :- xg) :- normalden(mu[j-1] :- xg)) :* -gama);
	}
	mes = (mes), (normalden(mu[ncat-1] :- xg) :* gama)
	
	//INCOMPLETE!!! if sumc(d) ne 0; .... incorporate dummies
	
	return(mes);
}

function op_me_raw(params, xbar, ncat, infcat, corresp,  | loop) {	 // loop as a parameter
	// infcat, corresp are not really needed. Added for compatibility
	k 		= cols(xbar)
	gama 	= params[1::k]
	mu 		= params[(k+1)::length(params)] 	// must be ordered!!!!!!!!!
	xg 		= xbar * gama
	
	mes		= normalden(mu[1] :- xg) :* -gama;	// ME for y=y0 with respect to all x' : k \times 1 vector
	
	for (j=2; j<=ncat-1; j++){
		mes = (mes) , ((normalden(mu[j] :- xg) :- normalden(mu[j-1] :- xg)) :* -gama);
	}
	mes = (mes), (normalden(mu[ncat-1] :- xg) :* gama)
	
	return(mes);
}

// derivative machine for ME in OP model
void _op_me_deriv(params, x, q, ncat, mevector) {
	me	= op_me(params', x, q, ncat)
	mevector	= rowshape(me,1)	// make a row vector
}


//


////////////////////////////////////////////////////////////////
// MIOPR = ZERO INFLATED ORDERED PROBIT MODEL with different X and Z variables, without correlation
////////////////////////////////////////////////////////////////

void _miopr_params(params, kx, kz, ncat, b, a, g, mu) {	
	b 	= params[1::kx]
	a	= params[kx+1]
	g	= params[(kx+2)::(kx+kz+1)]
	mu	= params[(kx+kz+2)::(kx+kz+1 + ncat - 1)]
	//mu	= punishSort(mu)
}


function MLmiopr(params, x, z, q, ncat, infcat, | loop) {
	kx 	= cols(x)
	kz 	= cols(z)
	n	= rows(x)
	
	_miopr_params(params, kx, kz, ncat, b = ., a = ., g = ., mu = .)
	
	xb	= x * b
	zg	= z * g
	p0	= normal(a :- xb)
	p1	= 1 :- p0
	
	prob = normal(mu'[J(n,1,1),] :- zg) , J(n,1,1)
	prob[ , 2::ncat] = prob[ , 2::ncat] - prob[ , 1::(ncat-1)]
	prob	= prob :* p1
	prob[ , infcat]	= prob[ , infcat] :+ p0
	// INCOMPLETE: check sum of probabilities
	
	if(loop == 1) {
		return(prob)
	} else if (loop == 2) {
		return((p0, prob[,infcat] - p0))
	} else if (loop == 3) {
		return ((p0, p1))
	}else{
		col_logl 	= rowsum(log(prob) :* q)
		return(col_logl)
	}
}

// d0 optimization machine for ZIOP model
void _miopr_optim(todo, params, x, z, q, ncat, infcat, coded, v, g, H) {
	kx 	= cols(x)
	kz 	= cols(z)
	n	= rows(x)
	_miopr_params(params', kx, kz, ncat, b = ., a = ., g = ., mu = .)
	if (coded == 1) {
		a  = decodeIncreasingSequence(a);
		mu = decodeIncreasingSequence(mu);
	}
	decoded_params = b \ a \ g \ mu 
	
	v = MLmiopr(decoded_params, x, z, q, ncat, infcat, 0)
	if(todo==1){
		// alas! gradient is not available so far
		//grad = Jacop(params', x, q, ncat)
		//g = colsum(grad);
	}
}


function miopr_me_raw(params, xzbar, ncat, infcat, corresp, |loop) {	 
	
	xbar = select(xzbar, corresp[1,]);
	zbar = select(xzbar, corresp[2,]);
	wherex = selectindex(corresp[1,]')
	wherez = selectindex(corresp[2,]')
	
	kx 	= cols(xbar)
	kz 	= cols(zbar)
	
	_miopr_params(params, kx, kz, ncat, b = ., a = ., g = ., mu = .)
	
	xb	= xbar * b
	zg	= zbar * g
	
	p_b	= 1 - normal(a - xb)
	p_g	= normal(mu' :- zg)
	p_g	= (p_g , 1) - (0, p_g)
	
	
	dp_b	= -normalden(a - xb) :* (-b)
	dp_g	= normalden(mu' :- zg)
	dp_g	= ((dp_g, 0) - (0, dp_g))[J(kz,1,1), ] :* (-g)
	
	if (loop == 2) {
		// marginal effects of probabilities of different zeroes
		me0x = (-dp_b, dp_b :* p_g[J(kx,1,1),infcat])
		me0z = (J(kz, 1, 0), p_b :* dp_g[,infcat])
		me0 = J(cols(corresp), 2, 0)
		me0[wherex, ] = me0[wherex, ] + me0x
		me0[wherez, ] = me0[wherez, ] + me0z
		return(me0)
	} else if (loop == 3) {
		mer = J(cols(corresp), 2, 0)
		mer[wherex, ] 	=  (-dp_b, dp_b)
		return(mer)
	}
	
	mex = dp_b :* p_g[J(kx,1,1), ]
	mex[, infcat] = mex[, infcat] :- dp_b
	mez = p_b :* dp_g
	
	
	mes = J(cols(corresp), cols(mex), 0)
	mes[wherex, ] = mes[wherex, ] + mex
	mes[wherez, ] = mes[wherez, ] + mez
	
	return(mes)
}








////////////////////////////////////////////////////////////////
// MIOPRc = ZERO INFLATED ORDERED PROBIT MODEL with different X and Z variables, with correlation
////////////////////////////////////////////////////////////////

void _mioprc_params(params, kx, kz, ncat, b, a, g, mu, ro) {	
	b 	= params[1::kx]
	a	= params[kx+1]
	g	= params[(kx+2)::(kx+kz+1)]
	mu	= params[(kx+kz+2)::(kx+kz+1 + ncat - 1)]
	ro	= params[ kx+kz+1 + ncat ]
}
function MLmioprc(params, x, z, q, ncat, infcat, | loop) {
	kx 	= cols(x)
	kz 	= cols(z)
	n	= rows(x)
	
	_mioprc_params(params, kx, kz, ncat, b = ., a = ., g = ., mu = ., ro = .)
	
	xb	= x * b
	zg	= z * g
	
	p0	= normal(a :- xb)
	p1	= 1 :- p0
	
	/* multiply 1'st and 3'rd arg by -1 because of ">=" inequality */
	prob = binormal(-(a:-xb), mu'[J(n,1,1),]:-zg, -ro), p1
	prob[,2::ncat] = prob[,2::ncat] - prob[,1::(ncat-1)]
	prob[ , infcat]	= prob[ , infcat] :+ p0
	
	if(loop == 1) {
		return(prob)
	} else if (loop == 2) {
		return((p0, prob[,infcat] - p0))
	} else if (loop == 3) {
		return ((p0, p1))
	}else{
		col_logl 	= rowsum(log(prob) :* q)
		return(col_logl)
	}
}
void _mioprc_optim(todo, params, x, z, q, ncat, infcat, coded, v, g, H) {
	kx 	= cols(x)
	kz 	= cols(z)
	n	= rows(x)
	_mioprc_params(params', kx, kz, ncat, b = ., a = ., g = ., mu = ., ro = .)
	if (coded == 1) {
		a  = decodeIncreasingSequence(a);
		mu = decodeIncreasingSequence(mu);
		ro = invlogit(ro) * 2 - 1
	}
	decoded_params = b \ a \ g \ mu \ ro
	
	v = MLmioprc(decoded_params, x, z, q, ncat, infcat, 0)
	if(todo==1){
		// alas! gradient is not available so far
	}
}
function mioprc_me_raw(params, xzbar, ncat, infcat, corresp, |loop) {
	xbar = select(xzbar, corresp[1,]);
	zbar = select(xzbar, corresp[2,]);
	wherex = selectindex(corresp[1,]')
	wherez = selectindex(corresp[2,]')
	
	kx 	= cols(xbar)
	kz 	= cols(zbar)
	
	_mioprc_params(params, kx, kz, ncat, b = ., a = ., g = ., mu = ., ro = .)
	
	xb	= xbar * b
	zg	= zbar * g
	
	// incomplete! Need to find that fucking derivatives, but I'll do it tomorrow (12 March)
	
	// Calculate partial derivatives of likelihood
	dp0 = normalden(a - xb) // normalden(xb-a) /* dP(0)/d(a-xb) */
	
	// start with derivates in non-inflated case
	/*
	p(y=5) = p(e1<xb-a, e2<zg-mu | ro)
	*/
	
	/*
	p2dx = ( dp0 :* normal(((mu' :- zg) :- ((xb :-a) * ro)) :/ sqrt(1 - ro^2)) , 0 )[J(kx,1,1),] :* b
	p2dx[,2::ncat] = p2dx[ ,2::ncat] :- p2dx[ , 1::(ncat-1)]
	*/
	p2dx = ( dp0 :* normal(((mu' :- zg) :- ((xb :-a) * ro)) :/ sqrt(1 - ro^2)) , dp0 )[J(kx,1,1),] :* b
	p2dx[,2::ncat] = p2dx[ ,2::ncat] :- p2dx[ , 1::(ncat-1)]
	
	p2dz = (normalden(mu' :- zg) :* normal(((xb:-a) :- ro * (mu' :- zg)) :/ sqrt(1 - ro^2)), 0)[J(kz,1,1),] :* -g
	p2dz[,2::ncat] = p2dz[ ,2::ncat] :- p2dz[ , 1::(ncat-1)]
	// add inflated case later
	
	
	// now calculate ME themselves
	
	if (loop == 2) {
		// marginal effects of probabilities of different zeroes
		me0x = ((dp0) :* (-b), p2dx[,infcat])
		me0z = (J(kz, 1, 0), p2dz[, infcat])
		
		me0 = J(cols(corresp), 2, 0)
		me0[wherex, ] 	= me0[wherex, ] + me0x
		me0[wherez, ] 	= me0[wherez, ] + me0z
		
		return(me0)
	} else if (loop == 3) {
		// marginal effects of zero / nonzero
		mer = J(cols(corresp), 2, 0)
		mer[wherex, ] 	=  (dp0, -dp0) :* -b[,(1,1,1)]
		return(mer)
	}
	
	mex	= p2dx
	mex[,infcat]	= mex[,infcat] + dp0 :* (-b)
	mez	= p2dz
	
	mes = J(cols(corresp), ncat, 0)
	mes[wherex, ] 	= mes[wherex, ] + mex
	mes[wherez, ] 	= mes[wherez, ] + mez
	
	return(mes)
	
}

////////////////////////////////////////////////////////////////
// CNOP model - without correlation
////////////////////////////////////////////////////////////////

void _cnop_params(params, kx, kzp, kzn, ncatp, ncatn, b, a, gp, mup, gn, mun) {
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp)]
	gn	= params[(kx+2+kzp+ncatp+1)::(kx+2+kzp+ncatp+kzn)]
	mun	= params[(kx+2+kzp+ncatp+kzn+1)::(kx+2+kzp+ncatp+kzn+ncatn)]
}

void _cnop_optim(todo, params, x, zp, zn, q, ncat, infcat, coded, v, g, H) {
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	_cnop_params(params', kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = .)
	if (coded == 1) {
		a  = decodeIncreasingSequence(a);
		mup = decodeIncreasingSequence(mup);
		mun = decodeIncreasingSequence(mun);
	} else {
		//"I AM THE LIKELIHOOD WRAPPER AND I DO NOT REPARAMETRIZE NOTHING!!!"
	}
	decoded_params = b \ a \ gp \ mup \ gn \ mun
	
	v = MLcnop(decoded_params, x, zp, zn, q, ncat, infcat, 0)
	//v = MLcnop(params', x, zp, zn, q, ncat, infcat, 0)
	if(todo==1) {
		// alas! gradient is not available so far
		grad = cnop_deriv(decoded_params', x, zp, zn, q, ncat, infcat, 1)
		g = grad;
	}
	if (sum(v:==.) > 0) {
		"Karramba! likelihood is null!"
	}
	//sum(v)
}

function MLcnop(params, x, zp, zn, q, ncat, infcat, | loop) {
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	_cnop_params(params, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = .)

	//a 	= punishSort(a)
	
	punishment = 0
	
	//a = punishSort(a, 10^-3, diff1=.)
	//mup = punishSort(mup, 10^-3, diff2=.)
	//mun = punishSort(mun, 10^-3, diff3=.)
	//punishment = punishment + diff1'*diff1 + diff2'*diff2 + diff3'*diff3
	//punishment = punishment * 100
	
	
	// INCOMPLETE: check dimensions
	xb = x*b
	zgp = zp*gp
	zgn = zn*gn
	
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	p2p = normal(mup'[J(n,1,1),] :- zgp) , J(n,1,1)
	p2p[ ,2::(ncatp+1)] = p2p[,2::(ncatp+1)] - p2p[,1::ncatp]
	
	p2n = normal(mun'[J(n,1,1),] :- zgn) , J(n,1,1)
	p2n[,2::(ncatn+1)] = p2n[,2::(ncatn+1)] - p2n[,1::ncatn]

	prob = J(n,ncat,0)
	prob[,infcat] = prob[,infcat]+p10
	prob[,1::(ncatn+1)] = prob[,1::(ncatn+1)] + p2n:*p1n
	prob[,infcat::ncat] = prob[,infcat::ncat] + p2p:*p1p
	
	// INCOMPLETE: check sum of probabilities
	
	if (loop == 1) {
		return(prob)
	} else if (loop == 2) {
		return((p10, p2n[,infcat] :* p1n, p2p[,1] :* p1p))
	} else if (loop == 3) {
		return ((p1n, p10, p1p))
	} else {
		col_logl 	= rowsum(log(prob) :* q)
		return(col_logl :- punishment)
	}
}


function cnop_deriv(params, x, zp, zn, q, ncat, infcat,| ofLogLik) {
	if(cols(ofLogLik) == 0) {
		ofLogLik = 0
	} // indicates, whether d(L) or d(ln L) should be returned
	
	
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a 	= punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp)]
	mup = punishSort(mup)
	gn	= params[(kx+2+kzp+ncatp+1)::(kx+2+kzp+ncatp+kzn)]
	mun	= params[(kx+2+kzp+ncatp+kzn+1)::(kx+2+kzp+ncatp+kzn+ncatn)]
	mun = punishSort(mun)
	
	// INCOMPLETE: check dimensions
	xb = x*b
	zgp = zp*gp
	zgn = zn*gn
	
	// probabilities of all events
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	p2p = normal(mup'[J(n,1,1),] :- zgp) , J(n,1,1)
	p2p[ ,2::(ncatp+1)] = p2p[,2::(ncatp+1)] - p2p[,1::ncatp]
	
	p2n = normal(mun'[J(n,1,1),] :- zgn) , J(n,1,1)
	p2n[,2::(ncatn+1)] = p2n[,2::(ncatn+1)] - p2n[,1::ncatn]

	prob = J(n,ncat,0)
	prob[,infcat] = prob[,infcat]+p10
	prob[,1::(ncatn+1)] = prob[,1::(ncatn+1)] + p2n:*p1n
	prob[,infcat::ncat] = prob[,infcat::ncat] + p2p:*p1p
	
	
	isn = rowsum(q[, 1::infcat])
	isp = rowsum(q[, infcat::ncat])
	is0 = q[,infcat]
	
	//  probabilities of actual events
	
	proba= rowsum(prob:* q)
	p2pa = rowsum(p2p :* q[,infcat::ncat])
	p2na = rowsum(p2n :* q[,1::infcat])
	p1pa = rowsum(p1p :* isp)
	p1na = rowsum(p1n :* isn)
	p10a = rowsum(p10 :* is0)
	
	// useful pdf's
	p2pd = normalden(mup'[J(n,1,1),] :- zgp) , J(n,1,0)
	p2pd[ ,2::(ncatp+1)] = p2pd[,2::(ncatp+1)] - p2pd[,1::ncatp]
	p2pda = rowsum(p2pd :* q[,infcat::ncat])
	
	p2nd = normalden(mun'[J(n,1,1),] :- zgn) , J(n,1,0)
	p2nd[,2::(ncatn+1)] = p2nd[,2::(ncatn+1)] - p2nd[,1::ncatn]
	p2nda = rowsum(p2nd :* q[,1::infcat])
	
	grad = J(n, kx+2+kzp+ncatp+kzn+ncatn, 0)
	
	
	// by beta
	for(j=1; j<=kx; j++){
		c = j + 0
		grad[,c]= grad[,c]	+ normalden(a[2]:-xb) :* x[,j] :* p2pa :* isp
		grad[,c]= grad[,c]	- normalden(a[2]:-xb) :* x[,j] 		 :* is0 
		grad[,c]= grad[,c]	+ normalden(a[1]:-xb) :* x[,j]		 :* is0 
		grad[,c]= grad[,c]	- normalden(a[1]:-xb) :* x[,j] :* p2na :* isn
	}
	// by alpha
	grad[,kx+1]	= grad[,kx+1] + ( 1) :* normalden(a[1]:-xb) :* p2na :* isn
	grad[,kx+1]	= grad[,kx+1] + (-1) :* normalden(a[1]:-xb)		  :* is0
	grad[,kx+2]	= grad[,kx+2] + (-1) :* normalden(a[2]:-xb) :* p2pa :* isp
	grad[,kx+2]	= grad[,kx+2] + ( 1) :* normalden(a[2]:-xb)		  :* is0
	// by gamma
	for(j=1; j<=kzp; j++){
		c = j+kx+2 
		grad[,c]= grad[,c] + p1pa :* (-zp[,j]) :* p2pda
	}
	
	for(j=1; j<=kzn; j++){
		c = j+kx+2+kzp+ncatp
		grad[,c]= grad[,c] + p1na :* (-zn[,j]) :* p2nda
	}
	// by mu
	c = 1+kx+2+kzp
	for(j=1; j<=ncatp; j++){
		c = j+kx+2+kzp 
		grad[,c]= grad[,c] + p1pa :* normalden(mup[j] :- zgp) :* (q[,infcat+j-1]-q[,infcat+j])
	}
	for(j=1; j<=ncatn; j++){
		c = j+kx+2+kzp+ncatp+kzn
		grad[,c]= grad[,c] + p1na:* normalden(mun[j] :- zgn) :* (q[,j]-q[,j+1])
	}
	if (ofLogLik){
		return (grad :/ proba)
	} else {
		return (grad)
	}
}


function cnop_me_raw(params, xzbar, ncat, infcat, corresp, |loop) {
	xbar 	= select(xzbar, corresp[1,]);
	zpbar 	= select(xzbar, corresp[2,]);
	znbar 	= select(xzbar, corresp[3,]);
	wherex 	= selectindex(corresp[1,]')
	wherezp = selectindex(corresp[2,]')
	wherezn = selectindex(corresp[3,]')
	
	kx	= cols(xbar)
	kzp	= cols(zpbar)
	kzn	= cols(znbar)
	
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a 	= punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp)]
	mup = punishSort(mup)
	gn	= params[(kx+2+kzp+ncatp+1)::(kx+2+kzp+ncatp+kzn)]
	mun	= params[(kx+2+kzp+ncatp+kzn+1)::(kx+2+kzp+ncatp+kzn+ncatn)]
	mun = punishSort(mun)
	
	xb  = xbar*b
	zgp = zpbar*gp
	zgn = znbar*gn
		
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	p2p = normal(mup' :- zgp) , 1
	p2p[ ,2::(ncatp+1)] = p2p[ ,2::(ncatp+1)] - p2p[,1::ncatp]
	
	p2n = normal(mun' :- zgn) , 1
	p2n[,2::(ncatn+1)] = p2n[,2::(ncatn+1)] - p2n[,1::ncatn]
	
	// now calculate ME themselves
	
	dp_bp	= -normalden(a[2] - xb) :* (-b)
	dp_gp	= normalden(mup' :- zgp)
	dp_gp	= ((dp_gp, 0) - (0, dp_gp))[J(kzp,1,1), ] :* (-gp)
	
	dp_bn	= normalden(a[1] - xb) :* (-b)
	dp_gn	= normalden(mun' :- zgn)
	dp_gn	= ((dp_gn, 0) - (0, dp_gn))[J(kzn,1,1), ] :* (-gn)
	
	if (loop == 2) {
		// marginal effects of probabilities of different zeroes
		me0x = (-dp_bp-dp_bn, dp_bn :* p2n[J(kx,1,1),infcat], dp_bp :* p2p[J(kx,1,1),1])
		me0zn = (J(kzn, 1, 0), p1n :* dp_gn[,infcat], J(kzn, 1, 0))
		me0zp = (J(kzp, 2, 0), p1p :* dp_gp[,1])
		me0 = J(cols(corresp), 3, 0)
		me0[wherex, ] = me0[wherex, ] + me0x
		me0[wherezp, ] = me0[wherezp, ] + me0zp
		me0[wherezn, ] = me0[wherezn, ] + me0zn
		return(me0)
	} else if (loop == 3) {
		mer = J(cols(corresp), 3, 0)
		mer[wherex, ] 	=  (dp_bn, -dp_bp-dp_bn, dp_bp)
		return(mer)
	}
	
	mex	 = J(kx, ncat, 0)
	mex[,1::infcat] = mex[,1::infcat] + p2n[J(kx,1,1),] :* dp_bn
	mex[,infcat] = mex[,infcat] - dp_bp - dp_bn
	mex[,infcat::ncat] = mex[,infcat::ncat] + p2p[J(kx,1,1),]  :* dp_bp
	
	mezp	= J(kzp, ncatn, 0), (p1p :* dp_gp)
	mezn	= (p1n :* dp_gn), J(kzn, ncatp, 0)
	
	
	mes = J(cols(corresp), ncat, 0)
	mes[wherex, ] 	= mes[wherex, ] + mex
	mes[wherezp, ] 	= mes[wherezp, ] + mezp
	mes[wherezn, ] 	= mes[wherezn, ] + mezn
	
	return(mes)
}

////////////////////////////////////////////////////////////////
// CNOPC model - with correlation
////////////////////////////////////////////////////////////////
void _cnopc_params(params, kx, kzp, kzn, ncatp, ncatn, b, a, gp, mup, gn, mun, rop, ron) {
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp)]
	gn	= params[(kx+2+kzp+ncatp+1)::(kx+2+kzp+ncatp+kzn)]
	mun	= params[(kx+2+kzp+ncatp+kzn+1)::(kx+2+kzp+ncatp+kzn+ncatn)]
	rop = params[(kx+2+kzp+ncatp+kzn+ncatn+1)]
	ron = params[(kx+2+kzp+ncatp+kzn+ncatn+2)]
}
void _cnopc_optim(todo, params, x, zp, zn, q, ncat, infcat, coded, lambda, v, g, H) {
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	_cnopc_params(params', kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = ., rop = ., ron = .)
	if (coded == 1) {
		a  = decodeIncreasingSequence(a);
		mup = decodeIncreasingSequence(mup);
		mun = decodeIncreasingSequence(mun);
		ron = invlogit(ron) * 2 - 1
		rop = invlogit(rop) * 2 - 1
	}
	decoded_params = b \ a \ gp \ mup \ gn \ mun \ rop \ ron
	
	v = MLcnopc(decoded_params, x, zp, zn, q, ncat, infcat, 0)
	
	/* 
	Add L2 regularization, making all params closer to 0 in optimum 
	Because v is vector of n observation, penalize each by lambda/n 
	*/
	//v = v :- lambda / n * sum(params :^ 2)
	
	if(todo==1){
		// alas! gradient is not available so far
	}
}


function MLcnopc(params, x, zp, zn, q, ncat, infcat, | loop) {
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	_cnopc_params(params, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = ., rop = ., ron = .)
	
	a = punishSort(a)
	mup = punishSort(mup)
	mun = punishSort(mun)
	
	corrmax = 0.9999
	punishment = 0;
	/*
	punishment = -1*(max((rop-1, 0, -rop-1))^2 + max((ron-1,0,-ron-1))^2)
	
	if(ron >= corrmax)	{ron = corrmax;	}
	if(ron <= -corrmax)	{ron = -corrmax;	}
	if(rop >= corrmax)	{rop = corrmax;	}
	if(rop <= -corrmax)	{rop = -corrmax;	}
	*/
	// INCOMPLETE: check dimensions
	xb = x*b
	zgp = zp*gp
	zgn = zn*gn
	
	// from here IT MUST BE CHANGED
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	
	pn = binormal(a[1]:-xb, mun'[J(n,1,1),] :- zgn, ron), p1n
	pn[,2::(ncatn+1)] = pn[ ,2::(ncatn+1)] :- pn[ , 1::ncatn]
	pp = binormal(-a[2]:+xb, mup'[J(n,1,1),] :- zgp, -rop), p1p
	pp[ ,2::(ncatp+1)] = pp[ ,2::(ncatp+1)] - pp[ ,1::ncatp]
	
	prob = J(n,ncat,0)
	prob[,infcat] = p10
	prob[,1::(ncatn+1)] = prob[,1::(ncatn+1)] + pn
	prob[,infcat::ncat] = prob[,infcat::ncat] + pp
	
	// INCOMPLETE: check sum of probabilities
	
	
	if (loop == 1) {
		return(prob)
	} else if (loop == 2) {
		return((p10, pn[,infcat], pp[,1]))
	} else if (loop == 3) {
		return ((p1n, p10, p1p))
	} else {
		col_logl 	= rowsum(log(prob) :* q) :+ punishment
		return(col_logl)
	}
}


function cnopc_me_raw(params, xzbar, ncat, infcat, corresp, |loop) {

	xbar 	= select(xzbar, corresp[1,]);
	zpbar 	= select(xzbar, corresp[2,]);
	znbar 	= select(xzbar, corresp[3,]);
	wherex 	= selectindex(corresp[1,]')
	wherezp = selectindex(corresp[2,]')
	wherezn = selectindex(corresp[3,]')
	
	kx	= cols(xbar)
	kzp	= cols(zpbar)
	kzn	= cols(znbar)
	
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp)]
	gn	= params[(kx+2+kzp+ncatp+1)::(kx+2+kzp+ncatp+kzn)]
	mun	= params[(kx+2+kzp+ncatp+kzn+1)::(kx+2+kzp+ncatp+kzn+ncatn)]
	rop = params[(kx+2+kzp+ncatp+kzn+ncatn+1)]
	ron = params[(kx+2+kzp+ncatp+kzn+ncatn+2)]

	corrmax = 0.9999
	if(ron >= corrmax)	{ron = corrmax;	}
	if(ron <= -corrmax)	{ron = -corrmax;	}
	if(rop >= corrmax)	{rop = corrmax;	}
	if(rop <= -corrmax)	{rop = -corrmax;	}
	
	xb = xbar*b
	zgp = zpbar*gp
	zgn = znbar*gn
	
	dp1p = normalden(a[1] - xb)
	dp1n = -normalden(a[2] - xb)
	
	// Calculate partial derivatives of likelihood
	
	pndx = ( dp1p :* normal(((mun' :- zgn) :- ((a[1] :-xb) * ron)) :/ sqrt(1 - ron^2)) , dp1p )[J(kx,1,1),] :* -b
	pndx[,2::(ncatn+1)] = pndx[ ,2::(ncatn+1)] :- pndx[ , 1::ncatn]
	pndz = (normalden(mun' :- zgn) :* normal(((a[1]:-xb) :- ron * (mun' :- zgn)) :/ sqrt(1 - ron^2)), 0)[J(kzn,1,1),] :* -gn
	pndz[,2::(ncatn+1)] = pndz[ ,2::(ncatn+1)] :- pndz[ , 1::ncatn]
	
	ppdx = ( dp1n:* normal(((mup' :- zgp) :+ ((-a[2]:+xb) * rop)) :/ sqrt(1 - rop^2)) , dp1n )[J(kx,1,1),] :* -b
	ppdx[ ,2::(ncatp+1)] = ppdx[ ,2::(ncatp+1)] - ppdx[ ,1::ncatp]
	ppdz = (normalden(mup' :- zgp) :* normal(((-a[2]:+xb) :+ rop * (mup' :- zgp)) :/ sqrt(1 - rop^2)), 0)[J(kzp,1,1),] :* -gp
	ppdz[ ,2::(ncatp+1)] = ppdz[ ,2::(ncatp+1)] - ppdz[ ,1::ncatp]

	// now calculate ME themselves
	
	if (loop == 2) {
		// marginal effects of probabilities of different zeroes
		me0x = ((-dp1p - dp1n) :* (-b), pndx[,infcat], ppdx[,1])
		me0zn = (J(kzn, 1, 0), pndz[, infcat], J(kzn, 1, 0))
		me0zp = (J(kzp, 2, 0), ppdz[,1])
		
		me0 = J(cols(corresp), 3, 0)
		me0[wherex, ] 	= me0[wherex, ] + me0x
		me0[wherezp, ] 	= me0[wherezp, ] + me0zp
		me0[wherezn, ] 	= me0[wherezn, ] + me0zn
		return(me0)
	} else if (loop == 3) {
		mer = J(cols(corresp), 3, 0)
		mer[wherex, ] 	=  (dp1n, -dp1p-dp1n, dp1p) :* -b[,(1,1,1)]
		return(mer)
	}
	
	
	mex	 = J(kx, ncat, 0)
	mex[,1::infcat] = mex[,1::infcat] + pndx
	mex[,infcat] = mex[,infcat] + (-dp1n - dp1p) :* (-b)
	mex[,infcat::ncat] = mex[,infcat::ncat] + ppdx
	
	mezp	= J(kzp, ncatn, 0), ppdz
	mezn	= pndz, J(kzp, ncatp, 0)
	
	mes = J(cols(corresp), ncat, 0)
	mes[wherex, ] 	= mes[wherex, ] + mex
	mes[wherezp, ] 	= mes[wherezp, ] + mezp
	mes[wherezn, ] 	= mes[wherezn, ] + mezn
	
	return(mes)
}


function cnopc_deriv(params, x, zp, zn, q, ncat, infcat) {
	// use the fact that d Binormal(x,y,ro)/d ro = exp(-1/2*(x^2+y^2-2*ro*x*y)/(1-ro^2))/(2*pi*sqrt(1-ro^2))
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a = punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp)]
	mup = punishSort(mup)
	gn	= params[(kx+2+kzp+ncatp+1)::(kx+2+kzp+ncatp+kzn)]
	mun	= params[(kx+2+kzp+ncatp+kzn+1)::(kx+2+kzp+ncatp+kzn+ncatn)]
	mun = punishSort(mun)
	rop = params[(kx+2+kzp+ncatp+kzn+ncatn+1)]
	ron = params[(kx+2+kzp+ncatp+kzn+ncatn+2)]
	
	corrmax = 0.9999
	
	if(ron >= corrmax)	{ron = 	corrmax;	}
	if(ron <= -corrmax)	{ron = -corrmax;	}
	if(rop >= corrmax)	{rop = 	corrmax;	}
	if(rop <= -corrmax)	{rop = -corrmax;	}
	
	// probabilities
	
	xb = xbar*b
	zgp = zpbar*gp
	zgn = znbar*gn
	
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	// matrix-form arguments to binormal CDF
	// width is 2-d, mup-d, mun-d 
	aa = (a[1]:-xb, a[2]:-xb)
	bbp= mup'[J(n,1,1),]:- zgp
	bbn= mun'[J(n,1,1),]:- zgn
	
	pn = binormal(aa[,1], bbp, ron), p1n
	pn[,2::(ncatn+1)] = pn[ ,2::(ncatn+1)] :- pn[ , 1::ncatn]
	pp = binormal(-aa[,2], bbn, -rop), p1p
	pp[ ,2::(ncatp+1)] = pp[ ,2::(ncatp+1)] - pp[ ,1::ncatp]
	
	
	prob = J(n,ncat,0)
	prob[,infcat] = p10
	prob[,1::(ncatn+1)] = prob[,1::(ncatn+1)] + pn
	prob[,infcat::ncat] = prob[,infcat::ncat] + pp
	
	proba= rowsum(prob:* q)
	
	isn = rowsum(q[, 1::infcat])
	isp = rowsum(q[, infcat::ncat])
	is0 = q[,infcat]
	
	// shifted left and right versions of label matrix (to mark neighbors)
	q_right = J(n, 1, 0), q[,2::ncat] 
	q_left = q[,1::(ncat-1)], J(n, 1, 0)
	
	// the gradient
	
	grad = J(n, kx+2+kzp+ncatp+kzn+ncatn+2, 0)
	
	// partial derivatives (also in matrix-form)
	// the fact is that dP(a, b, ro)/daa = fa*Fb
	// fa is some density, Fb is some probability, calculated below
	// but this 'fact' needs verification!
	fbbp	= normalden(bbp)
	fbbn 	= normalden(bbn)
	faa 	= normalden(aa)
	Fbbn 	= normal((bbn :- ron*a[1]) :/ sqrt(1 - ron^2)) 
	Fbbp 	= normal((bbp :- rop*a[2]) :/ sqrt(1 - rop^2)) 
	Faan 	= normal((a[1] :- ron*bbn) :/ sqrt(1 - ron^2)) 
	Faap 	= normal((-a[2] :+ rop*bbp):/ sqrt(1 - rop^2)) 
	
	
	
	/*
	for (i = 1; i <= n; i++) {
		dpda = 0
		dpdbp = 0
		dpdbn = 0
		dpdrp = 0
		dpdrn = 0
		
		
	}
	*/
	
	// INCOMPLETE!!!! all needs replacement
		// by beta
	for(j=1; j<=kx; j++){
		c = j + 0
		
		grad[,c]= grad[,c]	- faa[,2] :* x[,j] :* rowsum(((Fbbp, J(n,1,1)) - (J(n,1,0),Fbbp )) :* q[,infcat::ncatp]) 
		grad[,c]= grad[,c]	- normalden(a[2]:-xb) :* x[,j] 		 :* is0 
		
		grad[,c]= grad[,c]	+ faa[,1] :* x[,j] :* rowsum(((Fbbn, J(n,1,1)) - (J(n,1,0),Fbbn )) :* q[,1::infcat]) 
		grad[,c]= grad[,c]	+ normalden(a[1]:-xb) :* x[,j]			:* is0 
	}
	// by alpha
	grad[,kx+1]	= grad[,kx+1] + faa[,1] :* rowsum(q[,1::ncatn] 		:* ((Fbbn, J(n,1,1)) - (J(n,1,0),Fbbn)))
	grad[,kx+2]	= grad[,kx+2] - faa[,2] :* rowsum(q[,infcat::ncatp] :* ((Fbbp, J(n,1,1)) - (J(n,1,0),Fbbp))) 
	grad[,kx+1]	= grad[,kx+1] + (-1) :* normalden(a[1]:-xb)		  	:* is0	// need checking
	grad[,kx+2]	= grad[,kx+2] + ( 1) :* normalden(a[2]:-xb)		  	:* is0
	// by gamma
	for(j=1; j<=kzp; j++){
		c = j+kx+2 
		//grad[,c]= grad[,c] + p1pa :* (-zp[,j]) :* p2pda
	}
	
	for(j=1; j<=kzn; j++){
		c = j+kx+2+kzp+ncatp
		grad[,c]= grad[,c] + p1na :* (-zn[,j]) :* p2nda
	}
	// by mu
	c = 1+kx+2+kzp
	for(j=1; j<=ncatp; j++){
		c = j+kx+2+kzp 
		grad[,c]= grad[,c] + p1pa :* normalden(mup[j] :- zgp) :* (q[,infcat+j-1]-q[,infcat+j])
	}
	for(j=1; j<=ncatn; j++){
		c = j+kx+2+kzp+ncatp+kzn
		grad[,c]= grad[,c] + p1na:* normalden(mun[j] :- zgn) :* (q[,j]-q[,j+1])
	}
	// by rho
	grad[,kx+2+kzp+ncatp+kzn+ncatn+1] = 0
	grad[,kx+2+kzp+ncatp+kzn+ncatn+2] = 0
	
	if (ofLogLik){
		return (grad :/ proba)
	} else {
		return (grad)
	}
	
	
}





////////////////////////////////////////////////////////////////
// NOP model (without intersection) // add constraints!
////////////////////////////////////////////////////////////////
void _nop_params(params, kx, kzp, kzn, ncatp, ncatn, b, a, gp, mup, gn, mun) {
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp-1)]
	gn	= params[(kx+2+kzp+ncatp)::(kx+2+kzp+ncatp+kzn-1)]
	mun	= params[(kx+2+kzp+ncatp+kzn)::(kx+2+kzp+ncatp+kzn+ncatn-2)]
	
}
void _nop_optim(todo, params, x, zp, zn, q, ncat, infcat, coded, v, g, H){
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	ncatp = ncat - infcat
	ncatn = infcat - 1
	_nop_params(params', kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = .)
	if (coded == 1) {
		a  = decodeIncreasingSequence(a);
		mup = decodeIncreasingSequence(mup);
		mun = decodeIncreasingSequence(mun);
	} else {
		//"I AM THE LIKELIHOOD WRAPPER AND I DO NOT REPARAMETRIZE NOTHING!!!"
	}
	decoded_params = b \ a \ gp \ mup \ gn \ mun

	v = MLnop(decoded_params, x, zp, zn, q, ncat, infcat, 0)
	if(todo==1){
		// alas! gradient is not available so far
		grad = nop_deriv(decoded_params, x, zp, zn, q, ncat, infcat, 1)
		g = grad;
	}
}


function MLnop(params, x, zp, zn, q, ncat, infcat, | loop){
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a 	= punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp-1)]
	gn	= params[(kx+2+kzp+ncatp)::(kx+2+kzp+ncatp+kzn-1)]
	mun	= params[(kx+2+kzp+ncatp+kzn)::(kx+2+kzp+ncatp+kzn+ncatn-2)]
	
	punishment = 0
	mu = punishSort(mu, 10^-3, diff1=.)
	mup = punishSort(mup, 10^-3, diff2=.)
	mun = punishSort(mun, 10^-3, diff3=.)
	//punishment = punishment + diff1'*diff1 + diff2'*diff2 + diff3'*diff3
	punishment = punishment * 100
	
	
	// INCOMPLETE: check dimensions
	xb = x*b
	zgp = zp*gp
	zgn = zn*gn
	
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	p2p = normal(mup'[J(n,1,1),] :- zgp) , J(n,1,1)
	p2p[ ,2::(ncatp)] = p2p[,2::(ncatp)] - p2p[,1::(ncatp-1)]
	
	p2n = normal(mun'[J(n,1,1),] :- zgn) , J(n,1,1)
	p2n[,2::(ncatn)] = p2n[,2::(ncatn)] - p2n[,1::(ncatn-1)]

	prob = J(n,ncat,0)
	prob[,infcat] = p10
	prob[,1::ncatn] = p2n:*p1n
	prob[,(infcat+1)::ncat] = p2p:*p1p
	
	// INCOMPLETE: check sum of probabilities
	
	if (loop == 1) {
		return(prob)
	} else if (loop == 2) {
		// decomposition of zeros --- irrelevant
		//return((p10, p2n[,infcat] :* p1n, p2p[,1] :* p1p))
	} else if (loop == 3) {
		return ((p1n, p10, p1p))
	} else {
		col_logl 	= rowsum(log(prob) :* q)
		return(col_logl :- punishment)
	}
}

function nop_deriv(params, x, zp, zn, q, ncat, infcat,| ofLogLik) {
	if(cols(ofLogLik) == 0) {
		ofLogLik = 0
	} // indicates, whether d(L) or d(ln L) should be returned
	
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a 	= punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp-1)]
	gn	= params[(kx+2+kzp+ncatp)::(kx+2+kzp+ncatp+kzn-1)]
	mun	= params[(kx+2+kzp+ncatp+kzn)::(kx+2+kzp+ncatp+kzn+ncatn-2)]
	
	xb = x*b
	zgp = zp*gp
	zgn = zn*gn
	
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	p2p = normal(mup'[J(n,1,1),] :- zgp) , J(n,1,1)
	p2p[ ,2::(ncatp)] = p2p[,2::(ncatp)] - p2p[,1::(ncatp-1)]
	
	p2n = normal(mun'[J(n,1,1),] :- zgn) , J(n,1,1)
	p2n[,2::(ncatn)] = p2n[,2::(ncatn)] - p2n[,1::(ncatn-1)]

	prob = J(n,ncat,0)
	prob[,infcat] = p10
	prob[,1::ncatn] = p2n:*p1n
	prob[,(infcat+1)::ncat] = p2p:*p1p
	
	//
	
	isn = rowsum(q[, 1::(infcat-1)])
	isp = rowsum(q[, (infcat+1)::ncat])
	is0 = q[,infcat]
	
	//  probabilities of actual events
	
	proba= rowsum(prob:* q)
	p2pa = rowsum(p2p :* q[,(infcat+1)::ncat])
	p2na = rowsum(p2n :* q[,1::(infcat-1)])
	p1pa = rowsum(p1p :* isp)
	p1na = rowsum(p1n :* isn)
	p10a = rowsum(p10 :* is0)
	
	// useful pdf's
	p2pd = normalden(mup'[J(n,1,1),] :- zgp) , J(n,1,0)
	p2pd[ ,2::(ncatp)] = p2pd[,2::(ncatp)] - p2pd[,1::(ncatp-1)]
	p2pda = rowsum(p2pd :* q[,(infcat+1)::ncat])
	
	p2nd = normalden(mun'[J(n,1,1),] :- zgn) , J(n,1,0)
	p2nd[,2::(ncatn)] = p2nd[,2::(ncatn)] - p2nd[,1::(ncatn-1)]
	p2nda = rowsum(p2nd :* q[,1::(infcat-1)])
	
	grad = J(n, kx+2+kzp+ncatp-1+kzn+ncatn-1, 0)
	
	
	// by beta
	for(j=1; j<=kx; j++){
		c = j + 0
		grad[,c]= grad[,c]	+ normalden(a[2]:-xb) :* x[,j] :* p2pa :* isp
		grad[,c]= grad[,c]	- normalden(a[2]:-xb) :* x[,j] 		 :* is0 
		grad[,c]= grad[,c]	+ normalden(a[1]:-xb) :* x[,j]		 :* is0 
		grad[,c]= grad[,c]	- normalden(a[1]:-xb) :* x[,j] :* p2na :* isn
	}
	// by alpha
	grad[,kx+1]	= grad[,kx+1] + ( 1) :* normalden(a[1]:-xb) :* p2na :* isn
	grad[,kx+1]	= grad[,kx+1] + (-1) :* normalden(a[1]:-xb)		  :* is0
	grad[,kx+2]	= grad[,kx+2] + (-1) :* normalden(a[2]:-xb) :* p2pa :* isp
	grad[,kx+2]	= grad[,kx+2] + ( 1) :* normalden(a[2]:-xb)		  :* is0
	// by gamma
	for(j=1; j<=kzp; j++){
		c = j+kx+2 
		grad[,c]= grad[,c] + p1pa :* (-zp[,j]) :* p2pda
	}
	
	for(j=1; j<=kzn; j++){
		c = j+kx+2+kzp+ncatp-1
		grad[,c]= grad[,c] + p1na :* (-zn[,j]) :* p2nda
	}
	// by mu
	c = 1+kx+2+kzp
	for(j=1; j<=ncatp-1; j++){
		c = j+kx+2+kzp 
		grad[,c]= grad[,c] + p1pa :* normalden(mup[j] :- zgp) :* (q[,infcat+j-1]-q[,infcat+j])
	}
	for(j=1; j<=ncatn-1; j++){
		c = j+kx+2+kzp+ncatp-1+kzn
		grad[,c]= grad[,c] + p1na:* normalden(mun[j] :- zgn) :* (q[,j]-q[,j+1])
	}
	if (ofLogLik){
		return (grad :/ proba)
	} else {
		return (grad)
	}
}


function nop_me_raw(params, xzbar, ncat, infcat, corresp, |loop) {
	xbar 	= select(xzbar, corresp[1,]);
	zpbar 	= select(xzbar, corresp[2,]);
	znbar 	= select(xzbar, corresp[3,]);
	wherex 	= selectindex(corresp[1,]')
	wherezp = selectindex(corresp[2,]')
	wherezn = selectindex(corresp[3,]')
	
	kx	= cols(xbar)
	kzp	= cols(zpbar)
	kzn	= cols(znbar)
	
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a 	= punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp-1)]
	gn	= params[(kx+2+kzp+ncatp)::(kx+2+kzp+ncatp+kzn-1)]
	mun	= params[(kx+2+kzp+ncatp+kzn)::(kx+2+kzp+ncatp+kzn+ncatn-2)]
	
	xb  = xbar*b
	zgp = zpbar*gp
	zgn = znbar*gn
		
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	p2p = normal(mup'[1,] :- zgp) , 1
	p2p[ ,2::(ncatp)] = p2p[,2::(ncatp)] - p2p[,1::(ncatp-1)]
	
	p2n = normal(mun'[1,] :- zgn) , 1
	p2n[,2::(ncatn)] = p2n[,2::(ncatn)] - p2n[,1::(ncatn-1)]
	
	// now calculate ME themselves
	
	dp_bp	= -normalden(a[2] - xb) :* (-b)
	dp_gp	= normalden(mup' :- zgp)
	dp_gp	= ((dp_gp, 0) - (0, dp_gp))[J(kzp,1,1), ] :* (-gp)
	
	dp_bn	= normalden(a[1] - xb) :* (-b)
	dp_gn	= normalden(mun' :- zgn)
	dp_gn	= ((dp_gn, 0) - (0, dp_gn))[J(kzn,1,1), ] :* (-gn)
	
	if (loop == 2) {
		// marginal effects of probabilities of different zeroes --- IRRELEVANT
		me0x = (-dp_bp-dp_bn, dp_bn :* p2n[J(kx,1,1),infcat], dp_bp :* p2p[J(kx,1,1),1])
		me0zn = (J(kzn, 1, 0), p1n :* dp_gn[,infcat], J(kzn, 1, 0))
		me0zp = (J(kzp, 2, 0), p1p :* dp_gp[,1])
		me0 = J(cols(corresp), 3, 0)
		me0[wherex, ] = me0[wherex, ] + me0x
		me0[wherezp, ] = me0[wherezp, ] + me0zp
		me0[wherezn, ] = me0[wherezn, ] + me0zn
		return(me0)
	} else if (loop == 3) {
		mer = J(cols(corresp), 3, 0)
		mer[wherex, ] 	=  (dp_bn, -dp_bp-dp_bn, dp_bp)
		return(mer)
	}
	
	mex	 = J(kx, ncat, 0)
	mex[,1::(infcat-1)] = mex[,1::(infcat-1)] + p2n[J(kx,1,1),] :* dp_bn
	mex[,infcat] = mex[,infcat] - dp_bp - dp_bn
	mex[,(infcat+1)::ncat] = mex[,(infcat+1)::ncat] + p2p[J(kx,1,1),]  :* dp_bp
	
	mezp	= J(kzp, ncatn+1, 0), (p1p :* dp_gp)
	mezn	= (p1n :* dp_gn), J(kzn, ncatp+1, 0)
	
	
	mes = J(cols(corresp), ncat, 0)
	mes[wherex, ] 	= mes[wherex, ] + mex
	mes[wherezp, ] 	= mes[wherezp, ] + mezp
	mes[wherezn, ] 	= mes[wherezn, ] + mezn
	
	return(mes)
}
////////////////////////////////////////////////////////////////
// NOP-C model// add constraints!
////////////////////////////////////////////////////////////////
void _nopc_params(params, kx, kzp, kzn, ncatp, ncatn, b, a, gp, mup, gn, mun, rop, ron) {
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp-1)]
	gn	= params[(kx+2+kzp+ncatp)::(kx+2+kzp+ncatp+kzn-1)]
	mun	= params[(kx+2+kzp+ncatp+kzn)::(kx+2+kzp+ncatp+kzn+ncatn-2)]
	rop = params[(kx+2+kzp+ncatp+kzn+ncatn-1)]
	ron = params[(kx+2+kzp+ncatp+kzn+ncatn-0)]
}
void _nopc_optim(todo, params, x, zp, zn, q, ncat, infcat, coded, v, g, H){
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	ncatp = ncat - infcat
	ncatn = infcat - 1
	_nopc_params(params', kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = ., rop = ., ron = .)
	if (coded == 1) {
		a  = decodeIncreasingSequence(a);
		mup = decodeIncreasingSequence(mup);
		mun = decodeIncreasingSequence(mun);
		ron = invlogit(ron) * 2 - 1;
		rop = invlogit(rop) * 2 - 1;
	} else {
		//"I AM THE LIKELIHOOD WRAPPER AND I DO NOT REPARAMETRIZE NOTHING!!!"
	}
	decoded_params = b \ a \ gp \ mup \ gn \ mun \ rop \ ron
	v = MLnopc(decoded_params, x, zp, zn, q, ncat, infcat, 0)
	if(todo==1){
		// alas! gradient is not available so far
		//grad = nopc_deriv(decoded_params, x, zp, zn, q, ncat, infcat, 1)
		//g = grad;
	}
}
function MLnopc(params, x, zp, zn, q, ncat, infcat, | loop){
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a 	= punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp-1)]
	gn	= params[(kx+2+kzp+ncatp)::(kx+2+kzp+ncatp+kzn-1)]
	mun	= params[(kx+2+kzp+ncatp+kzn)::(kx+2+kzp+ncatp+kzn+ncatn-2)]
	rop = params[(kx+2+kzp+ncatp+kzn+ncatn-1)]
	ron = params[(kx+2+kzp+ncatp+kzn+ncatn-0)]
	
	punishment = 0
	mu = punishSort(mu, 10^-3, diff1=.)
	mup = punishSort(mup, 10^-3, diff2=.)
	mun = punishSort(mun, 10^-3, diff3=.)
	//punishment = punishment + diff1'*diff1 + diff2'*diff2 + diff3'*diff3
	
	corrmax = 0.999;
	punishment = punishment + (max((rop-1, 0, -rop-1))^2 + max((ron-1,0,-ron-1))^2)
	
	if(ron >= 1)	{ron = corrmax;	}
	if(ron <= -1)	{ron = -corrmax;	}
	if(rop >= 1)	{rop = corrmax;	}
	if(rop <= -1)	{rop = -corrmax;	}
	
	punishment = punishment * 100000;
	
	
	// INCOMPLETE: check dimensions
	xb = x*b
	zgp = zp*gp
	zgn = zn*gn
	
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	
	pn = binormal(a[1]:-xb, mun'[J(n,1,1),] :- zgn, ron), p1n
	pn[,2::(ncatn)] = pn[ ,2::(ncatn)] :- pn[ , 1::(ncatn-1)]
	pp = binormal(-a[2]:+xb, mup'[J(n,1,1),] :- zgp, -rop), p1p
	pp[ ,2::(ncatp)] = pp[ ,2::(ncatp)] - pp[ ,1::(ncatp-1)]
	
	
	prob = J(n,ncat,0)
	prob[,infcat] = p10
	prob[,1::(ncatn)] = prob[,1::(ncatn)] + pn
	prob[,(infcat+1)::ncat] = prob[,(infcat+1)::ncat] + pp
	
	
	// INCOMPLETE: check sum of probabilities
	
	if (loop == 1) {
		return(prob)
	} else if (loop == 2) {
		// decomposition of zeros --- irrelevant
		//return((p10, p2n[,infcat] :* p1n, p2p[,1] :* p1p))
	} else if (loop == 3) {
		return ((p1n, p10, p1p))
	} else {
		col_logl 	= rowsum(log(prob) :* q)
		return(col_logl :- punishment)
	}
}

function nopc_deriv(params, x, zp, zn, q, ncat, infcat,| ofLogLik) {
	n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a 	= punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp-1)]
	gn	= params[(kx+2+kzp+ncatp)::(kx+2+kzp+ncatp+kzn-1)]
	mun	= params[(kx+2+kzp+ncatp+kzn)::(kx+2+kzp+ncatp+kzn+ncatn-2)]
	rop = params[(kx+2+kzp+ncatp+kzn+ncatn-1)]
	ron = params[(kx+2+kzp+ncatp+kzn+ncatn-0)]
	
	punishment = 0
	mu = punishSort(mu, 10^-3, diff1=.)
	mup = punishSort(mup, 10^-3, diff2=.)
	mun = punishSort(mun, 10^-3, diff3=.)
	//punishment = punishment + diff1'*diff1 + diff2'*diff2 + diff3'*diff3
	
	corrmax = 0.999;
	punishment = punishment + (max((rop-1, 0, -rop-1))^2 + max((ron-1,0,-ron-1))^2)
	
	if(ron >= 1)	{ron = corrmax;	}
	if(ron <= -1)	{ron = -corrmax;	}
	if(rop >= 1)	{rop = corrmax;	}
	if(rop <= -1)	{rop = -corrmax;	}
	
	punishment = punishment * 100000;
	
	
	// INCOMPLETE: check dimensions
	xb = x*b
	zgp = zp*gp
	zgn = zn*gn
	
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n
	
	
	pn = binormal(a[1]:-xb, mun'[J(n,1,1),] :- zgn, ron), p1n
	pn[,2::(ncatn)] = pn[ ,2::(ncatn)] :- pn[ , 1::(ncatn-1)]
	pp = binormal(-a[2]:+xb, mup'[J(n,1,1),] :- zgp, -rop), p1p
	pp[ ,2::(ncatp)] = pp[ ,2::(ncatp)] - pp[ ,1::(ncatp-1)]
	
	
	prob = J(n,ncat,0)
	prob[,infcat] = p10
	prob[,1::(ncatn)] = prob[,1::(ncatn)] + pn
	prob[,(infcat+1)::ncat] = prob[,(infcat+1)::ncat] + pp

	
	
	// Beginning from here, rewrite everything for the correlated version! 
	
	
	isn = rowsum(q[, 1::(infcat-1)])
	isp = rowsum(q[, (infcat+1)::ncat])
	is0 = q[,infcat]
	
	//  probabilities of actual events
	
	proba= rowsum(prob:* q)
	p2pa = rowsum(p2p :* q[,(infcat+1)::ncat])
	p2na = rowsum(p2n :* q[,1::(infcat-1)])
	p1pa = rowsum(p1p :* isp)
	p1na = rowsum(p1n :* isn)
	p10a = rowsum(p10 :* is0)
	
	// useful pdf's
	p2pd = normalden(mup'[J(n,1,1),] :- zgp) , J(n,1,0)
	p2pd[ ,2::(ncatp)] = p2pd[,2::(ncatp)] - p2pd[,1::(ncatp-1)]
	p2pda = rowsum(p2pd :* q[,(infcat+1)::ncat])
	
	p2nd = normalden(mun'[J(n,1,1),] :- zgn) , J(n,1,0)
	p2nd[,2::(ncatn)] = p2nd[,2::(ncatn)] - p2nd[,1::(ncatn-1)]
	p2nda = rowsum(p2nd :* q[,1::(infcat-1)])
	
	grad = J(n, kx+2+kzp+ncatp-1+kzn+ncatn-1, 0)
	
	
	// by beta
	for(j=1; j<=kx; j++){
		c = j + 0
		grad[,c]= grad[,c]	+ normalden(a[2]:-xb) :* x[,j] :* p2pa :* isp
		grad[,c]= grad[,c]	- normalden(a[2]:-xb) :* x[,j] 		 :* is0 
		grad[,c]= grad[,c]	+ normalden(a[1]:-xb) :* x[,j]		 :* is0 
		grad[,c]= grad[,c]	- normalden(a[1]:-xb) :* x[,j] :* p2na :* isn
	}
	// by alpha
	grad[,kx+1]	= grad[,kx+1] + ( 1) :* normalden(a[1]:-xb) :* p2na :* isn
	grad[,kx+1]	= grad[,kx+1] + (-1) :* normalden(a[1]:-xb)		  :* is0
	grad[,kx+2]	= grad[,kx+2] + (-1) :* normalden(a[2]:-xb) :* p2pa :* isp
	grad[,kx+2]	= grad[,kx+2] + ( 1) :* normalden(a[2]:-xb)		  :* is0
	// by gamma
	for(j=1; j<=kzp; j++){
		c = j+kx+2 
		grad[,c]= grad[,c] + p1pa :* (-zp[,j]) :* p2pda
	}
	
	for(j=1; j<=kzn; j++){
		c = j+kx+2+kzp+ncatp-1
		grad[,c]= grad[,c] + p1na :* (-zn[,j]) :* p2nda
	}
	// by mu
	c = 1+kx+2+kzp
	for(j=1; j<=ncatp-1; j++){
		c = j+kx+2+kzp 
		grad[,c]= grad[,c] + p1pa :* normalden(mup[j] :- zgp) :* (q[,infcat+j-1]-q[,infcat+j])
	}
	for(j=1; j<=ncatn-1; j++){
		c = j+kx+2+kzp+ncatp-1+kzn
		grad[,c]= grad[,c] + p1na:* normalden(mun[j] :- zgn) :* (q[,j]-q[,j+1])
	}
	
	// by rho
	
	
	
	
	// return
	if (ofLogLik){
		return (grad :/ proba)
	} else {
		return (grad)
	}
}


function nopc_me_raw(params, xzbar, ncat, infcat, corresp, |loop) {

	xbar 	= select(xzbar, corresp[1,]);
	zpbar 	= select(xzbar, corresp[2,]);
	znbar 	= select(xzbar, corresp[3,]);
	wherex 	= selectindex(corresp[1,]')
	wherezp = selectindex(corresp[2,]')
	wherezn = selectindex(corresp[3,]')
	
	kx	= cols(xbar)
	kzp	= cols(zpbar)
	kzn	= cols(znbar)
	
	// incomplete: assumption
	ncatp = ncat - infcat
	ncatn = infcat - 1
	
	b 	= params[1::kx]
	a	= params[(kx+1)::(kx+2)]
	a 	= punishSort(a)
	gp	= params[(kx+2+1)::(kx+2+kzp)]
	mup	= params[(kx+2+kzp+1)::(kx+2+kzp+ncatp-1)]
	gn	= params[(kx+2+kzp+ncatp)::(kx+2+kzp+ncatp+kzn-1)]
	mun	= params[(kx+2+kzp+ncatp+kzn)::(kx+2+kzp+ncatp+kzn+ncatn-2)]
	rop = params[(kx+2+kzp+ncatp+kzn+ncatn-1)]
	ron = params[(kx+2+kzp+ncatp+kzn+ncatn-0)]
	

	mu = punishSort(mu, 10^-3, diff1=.)
	mup = punishSort(mup, 10^-3, diff2=.)
	mun = punishSort(mun, 10^-3, diff3=.)
	//punishment = punishment + diff1'*diff1 + diff2'*diff2 + diff3'*diff3
	corrmax = 0.999;
	if(ron >= 1)	{ron = corrmax;	}
	if(ron <= -1)	{ron = -corrmax;	}
	if(rop >= 1)	{rop = corrmax;	}
	if(rop <= -1)	{rop = -corrmax;	}
	
	xb = xbar*b
	zgp = zpbar*gp
	zgn = znbar*gn
	
	dp1p = normalden(a[1] - xb)
	dp1n = -normalden(a[2] - xb)
	
	// Calculate partial derivatives of likelihood
	
	pndx = ( dp1p :* normal(((mun' :- zgn) :- ((a[1] :-xb) * ron)) :/ sqrt(1 - ron^2)) , dp1p )[J(kx,1,1),] :* -b
	pndx[,2::(ncatn)] = pndx[ ,2::(ncatn)] :- pndx[ , 1::(ncatn-1)]
	pndz = (normalden(mun' :- zgn) :* normal(((a[1]:-xb) :- ron * (mun' :- zgn)) :/ sqrt(1 - ron^2)), 0)[J(kzn,1,1),] :* -gn
	pndz[,2::(ncatn)] = pndz[ ,2::(ncatn)] :- pndz[ , 1::(ncatn-1)]
	
	ppdx = ( dp1n:* normal(((mup' :- zgp) :+ ((-a[2]:+xb) * rop)) :/ sqrt(1 - rop^2)) , dp1n )[J(kx,1,1),] :* -b
	ppdx[ ,2::(ncatp)] = ppdx[ ,2::(ncatp)] - ppdx[ ,1::(ncatp-1)]
	ppdz = (normalden(mup' :- zgp) :* normal(((-a[2]:+xb) :+ rop * (mup' :- zgp)) :/ sqrt(1 - rop^2)), 0)[J(kzp,1,1),] :* -gp
	ppdz[ ,2::(ncatp)] = ppdz[ ,2::(ncatp)] - ppdz[ ,1::(ncatp-1)]

	// now calculate ME themselves
	
	if (loop == 2) {
		// marginal effects of probabilities of different zeroes
		"Probabilities of diffierent zeroes are not what you need with NOPC! Zeroes are just zeroes with it."
		return(.)
	} else if (loop == 3) {
		mer = J(cols(corresp), 3, 0)
		mer[wherex, ] 	=  (dp1n, -dp1p-dp1n, dp1p) :* -b[,(1,1,1)]
		return(mer)
	}
	
	
	mex	 = J(kx, ncat, 0)
	mex[,1::(infcat-1)] = mex[,1::(infcat-1)] + pndx
	mex[,infcat] = mex[,infcat] + (-dp1n - dp1p) :* (-b)
	mex[,(infcat+1)::ncat] = mex[,(infcat+1)::ncat] + ppdx
	
	mezp	= J(kzp, ncatn + 1, 0), ppdz
	mezn	= pndz, J(kzp, ncatp + 1, 0)
	
	mes = J(cols(corresp), ncat, 0)
	mes[wherex, ] 	= mes[wherex, ] + mex
	mes[wherezp, ] 	= mes[wherezp, ] + mezp
	mes[wherezn, ] 	= mes[wherezn, ] + mezn
	
	return(mes)
}


end

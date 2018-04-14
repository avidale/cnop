version 12

mata: mata clear

// define CNOPishModel struct
run C:\Users\Dave\YandexDisk\hsework\gauss-mata\sandbox\CNOPishModel_definition.ado
// define all functions for optimization and ME
run C:\Users\Dave\YandexDisk\hsework\gauss-mata\sandbox\gradients.ado


mata

x = ((1::5)\(5::1)) ,(1::10)
y = (1,1,2,1,2,2,2,1,2,2)'


n	= rows(x)
k	= cols(x)

allcat = uniqrows(y)
ncat = rows(allcat)

q = J(n,ncat,0)

for(i=1;i<=ncat; i++){
	q[.,i]=(y :== allcat[i])
}


start_mu	= invnormal(runningsum(mean(q))[1::ncat-1])';
start_b		= invsym(x'*x)*x'*y;
start_param	= start_b \ start_mu 

S = optimize_init()
optimize_init_argument(S, 1, x)
optimize_init_argument(S, 2, q)
optimize_init_argument(S, 3, ncat)
optimize_init_evaluator(S, &_op_optim())
optimize_init_evaluatortype(S, "d1")
optimize_init_params(S, (start_param'))
params = optimize(S)'
params

maxLik	= optimize_result_value(S)
grad 	= optimize_result_gradient(S)
covMat	= optimize_result_V(S)
retCode	= optimize_result_returncode(S)
//covMat_rob	= optimize_result_V_robust(S)		// in contrast with the source, here I get it automatically

g	= Jacop(params, x, q, ncat) // gradient for all observations

/*
ss	= J(rows(params), rows(params), 0)

for (i=1; i<=maxc(who);i++){
	sel = select(g, who == i)
    ss  = ss + sel' * sel; 
}

covMat_rob =covMat * ss * covMat;  /* esimated robust asymptotic variance */
*/
covMat_rob = covMat // clarify that issue


se		= sqrt(diagonal(covMat))
tstat	= abs(params :/ se)

se_rob		= sqrt(diagonal(covMat))
tstat_rob	= abs(params :/ se_rob)

AIC		= -2 * maxLik + 2 * rows(params) 
BIC		= -2 * maxLik + ln(n) * rows(params)
CAIC	= -2 * maxLik + (1 + ln(n)) * rows(params)
AICc	= AIC + 2 * rows(params) * (rows(params) + 1) / (n - rows(params) - 1)
HQIC	= -2 * maxLik + 2*rows(params)*ln(ln(n))

pred_prob	= MLop(params, x, q, ncat, 1)
// predicted option??

me		= op_me(params, x, q, ncat)
mese	= J(k* ncat,1, 0) // to be reshaped
mese_rob	= J(k* ncat,1, 0) // to be reshaped

// delta method!!!


// I want to find derivative of all ME's with respect to all parameters
D = deriv_init()
deriv_init_evaluator(D, &_op_me_deriv())
deriv_init_evaluatortype(D, "t")
deriv_init_argument(D, 1, x)
deriv_init_argument(D, 2, q)
deriv_init_argument(D, 3, ncat)
deriv_init_params(D, params') // a row vector
dydx = deriv(D, 1)
// in columns - parameters
// in rows - elements of rowshape(me[,1::cols(x)],1)

for(i = 1; i<=rows(dydx); i++){
	mese[i]	= dydx[i,] * covMat * dydx[i,]'
	mese_rob[i]	= dydx[i,] * covMat_rob * dydx[i,]'
}
mese	= colshape(mese, ncat)
mese_rob	= colshape(mese_rob, ncat)


met		= abs(me) :/ sqrt(mese)
met_rob	= abs(me) :/ sqrt(mese_rob)

met

/*
gama = params[1::2]
mu = params[3]
xg 		= x * gama */


end

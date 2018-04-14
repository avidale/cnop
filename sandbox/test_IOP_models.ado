version 12

mata: mata clear

// subroutines from .ado files will be saved in MATA librarie upon completion

// define CNOPModel struct
run C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\CNOPishModel_definition.ado
// define all auxiliary functions for optimization and ME
run C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\gradients.ado
// define functions for estimation
run C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\inflatedOP_estimation_routines.ado
// define functions for generating data
run C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\generate_CNOP.ado


mata

// generate ZIOP data
/*
n = 1000;	k = 4;	ncat = 5;	infcat = 3;
beta	= ( 2 \ 0 \ 0 \ 4);	a = -2
gamma	= ( 3 \ 2 \ 1 \ 0);	mu	= (-2 \ -1 \  1 \  2)
x = y = q = z = 0;
genZIOP(x, y, q, n, k, ncat, infcat, beta, a, gamma, mu)


// instance of CNOPModel class
OP_model = estimateOP(y, x)
OP_model.params, OP_model.se //, OP_model.se_rob
// OP_model.logLik
// colsum(OP_model.probabilities)/rows(x)
//OP_model.met
*/

/*
MIOP_model = estimateMIOP(y, x, 3, 0)
(beta \ a \gamma \mu), MIOP_model.params, MIOP_model.se, MIOP_model.t

MIOPR_model = estimateMIOPR(y, x, x, 3, 0)
(beta \ a \gamma \mu), MIOPR_model.params, MIOPR_model.se, MIOPR_model.t

MIOPR2_model = estimateMIOPR(y, x[.,(1,4)], x, 3, 0)
(beta[(1,4)] \ a \gamma \mu), MIOPR2_model.params, MIOPR2_model.se, MIOPR2_model.t

MIOPR3_model = estimateMIOPR(y, x[.,(1,4)], x[.,1::3], 3, 0)
(beta[(1,4)] \ a \gamma[1::3] \mu), MIOPR3_model.params, MIOPR3_model.se, MIOPR3_model.t
*/

// Vuongtest(OP_model, MIOPR_model, q)


/*
n = 1000;	k = 5;	ncat = 7;	infcat = 4;
beta	= ( 2 \ 0 \ 0 \ 4\ 1);	a = (-2 \ 1 );
gamma	= ( 3 \ 2 \ 1 \ 0 \-1);	mu	= (-3 \ 0 \  3)
x = y = q = 0;
genMIOPR(x, y, q, n, k, ncat, infcat, beta, a, gamma, mu, gamma, mu)
colsum(q)

params = beta\a\gamma\mu\gamma\mu
//prob = MLcnopc(params, x, x, x, q, ncat, infcat, 1)
//params = beta\a\gamma\mu\gamma\mu\0.3\-0.3
//prob = MLmioprc(params, x, x, x, q, ncat, infcat, 1)

CNOP_model = estimateCNOP(y, x, x, x, 4, 0)
params, CNOP_model.params, CNOP_model.t

//MIOPRC_model = estimateMIOPRC(y, x, x, x, 4)

//(params\0\0) , MIOPRC_model.params,  MIOPRC_model.t
*/

n = 1000; k = 5; ncat = 7; infcat = 4;
beta = (2 \ 3 \ 0 \ 0 \ 0); a = (-0.4 \ 1.2);
gammap = (0 \ 1 \ 2 \ 0 \ 0.1); mup = (-5 \ -2 \ 3 );
gamman = (0 \ 0 \ 3 \ 5 \ 0.2); mun = (-4 \ -1 \ 3);

/*
beta = (2 \ 3 \ 0 \ 0 \ 0); a = (-0.4 \ 1.2);
gammap = (0 \ 0 \ 2 \ 0 \ 0); mup = (-5 \ -2 \ 3 );
gamman = (0 \ 0 \ 0 \ 5 \ 0.2); mun = (-4 \ -1 \ 3);
*/

x=y=q=0
genCNOP(x, y, q, n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun)
colsum(q)
params = beta\a\gammap\mup\gamman\mun

CNOP_model = estimateCNOP(y, x, x, x, 4)
params, CNOP_model.params, CNOP_model.se


//test grad
/*
grad = cnop_deriv(params, x, x, x, q, ncat, infcat, 1)
lik0 = (MLcnop(params, x,x,x,q, ncat, infcat))
delta = 0.00001
gradNum = J(n, length(params),0)
for(i = 1; i <= length(params); i++) {
	params2 = params
	params2[i] = params[i] + delta
	lik1 = (MLcnop(params2, x,x,x,q, ncat, infcat))
	gradNum[,i] = (lik1 - lik0) :/ delta
}
// (grad - gradNum)[1::5,]
colsum(gradNum) - colsum(grad)
*/

// test CNOPC
/*
//CNOPC_model = estimateCNOPC(y, x, x, x, 4, 0)
//(params \ 0 \ 0), CNOPC_model.params, CNOPC_model.t

rp = 0.5; rn = -0.5;
genCNOPC(x, y, q, n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun, rn, rp)
params = params\rn\rp
CNOPC_model = estimateCNOPC(y, x, x, x, 4, 0)
params, CNOPC_model.params, CNOPC_model.se
*/

/*
// testing raw MEs
corresp = J(3,5,1)
xzbar = J(1,5,0)

mana = cnopc_me_raw(params, xzbar, 7, 4, corresp)
mana
delta = 0.00001
mnum = J(5,7,0)
for(i = 1; i <= 5; i++) {
	xzb = xzbar 
	xzb[i] = xzbar[i] + delta
	mnum[i,] = (MLcnopc(params, xzb, xzb, xzb, J(1,7,1/7) , 7, 4,1 ) - MLcnopc(params, xzbar, xzbar, xzbar, J(1,7,1/7) , 7, 4,1)) :/ delta
}
(mnum - mana) :/ mnum
*/

end

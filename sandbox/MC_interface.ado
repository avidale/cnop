version 12

// include file that in turn includes all needed functions to estimate models
//run C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\stata_wrappers.ado

cd "C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox"

mata

function setcorresp(class CNOPModel scalar model, corresp) {
	// workaround MATA not recognizning class of models
	// reason is that MATA allows full variable declarations only inside functions Oo
	model.corresp = corresp
}

// functions which work with CSV (without header)
function writeMatrix(filename, data) {
	unlink(filename)
	fh = fopen(filename, "w")
	for (i = 1; i <= rows(data); i++) {
		fput(fh, invtokens(strofreal(data[i,]), ","))
	}
	fclose(fh)
}

function readMatrix(filename) {
	fh = fopen(filename, "r")
	tok = tokeninit(",")
	if ((line=fget(fh))!=J(0,0,"")) {
		tokenset(tok, line)
		result = strtoreal(tokengetall(tok))
	}
	while ((line=fget(fh))!=J(0,0,"")) {
		tokenset(tok, line)
		result = result \ strtoreal(tokengetall(tok))
	}
	fclose(fh)
	return(result)
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



x = readMatrix("x_2000.csv")

void doSimulations(filesuffix, No, rep, meNo, mod, dgp, overlap, allvars, cv, quiet, seed, me_detailed_idx, converged_only, precise_only, robust, who, b, bse, mese, mese0, probse, probse0, par_true, me_true, me0_true, p_true, p0_true, label, est_time, converged, error_codes, precise) {

	/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       MANUAL INPUT      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
	{
	//rseed(seed)	// random state
	
	//No          = 1000;    /* choose sample size */
	//rep         = 10;   /* number of generated samples */
	//meNo        = 10;	// NUMBER OF OBS TO CALCULATE ME 
	// dgp: 1 - op, 2 - nop, 3 - tiop
	//dgp = "CNOP"
	// overlap: 1 - no, 2 - partial, 3 - complete
	//overlap = 2

	cJ     = 5;      /* enter the number of outcome categories */
	inf    = 3;      /* indicate which category of y (1,2,3,...) is to be inflated */
	xbar   = 1;      /* enter values of covariates (observations in rows) to compute MEs; set to 0 for zeros, 1 for population means  */
	dummiesVector  = (0, 0, 1); 
					 /* indicators of dummy variables for ME estimation */
	bound  = 1000;   /* set bounds for parameters in cml procedure */
	
	// cv critical values for confidence intervals 
	
	//quiet = 0 	// quiet estimation of models
	
	// converged_only - calculate summary statisitcs only for converged models

	/*                       load or generate artificial covariates                     */
	x1   = rnormal(No,1,2,1);
	x2   = rnormal(No,1,0,1);
	x3   = runiform(No,1);
	x3   = (x3 :> 0.70) :- (x3 :< 0.30);
	x4   = rnormal(No,1,0,1);
	/* parameters are calibrated to produce in case of 5 outcome categories approximately:
	y: -2: 7%, -1: 14%, 0: 58%, 1: 14%, 2: 7% observations
	Decomposition of zeros:   r=-1: 1/3, r=0: 1/3, r=1: 1/3 
	*/
	}

	if (dgp=="CNOP" || dgp=="CNOPC"){ /*  TIOP */
		if (overlap == 1) {  //@ no overlap @
			xv   = 1 \ 0 \ 0 ;   
			znv  = 0 \ 1 \ 0 ;
			zpv  = 0 \ 0 \ 1 ;

			beta    = (.6 ,.4, .8)';
			mu      = (0.95, 1.45)';
			betan   = (.2, .8, .3 )';
			betap   = (.4 ,.3, .9 )';
			mun     = (-1.22, 0.03)';
			mup     = (-.03, 1.18)';
		} else if (overlap == 2) { //@ partial overlap @
			xv	= 1 \ 1 \ 0 ;   
			znv     = 1 \ 0 \ 1 ;
			zpv     = 0 \ 1 \ 1 ;

			beta    = (.6, .4, .8)';
			mu      = (.9, 1.5)';
			betan   = (.2, .8, .3)';
			betap   = (.4, .3, .9)';
			mun     = (-.67, .36)';
			mup     = (.02, 1.28)';
		} else if (overlap == 3){  //@ complete overlap @
			xv   = 1 \ 1 \ 1 ;   
			znv  = 1 \ 1 \ 1 ;
			zpv  = 1 \ 1 \ 1 ;

			beta    = (.6, .4, .8 )';
			mu      = (0.85, 1.55)';
			betan   = (.2, .8, .3 )';
			betap   = (.4, .3, .9 )';
			mun     = (-1.2, 0.07)';
			mup     = (1.28, 2.5)';
			
		}
		
		beta        = select(beta,  xv);
		betan       = select(betan, znv);
		betap       = select(betap, zpv);
		
		if (dgp=="CNOPC") {
			rop = 0.6;
			ron = -0.1;
		} else {
			rop = 0;
			ron = 0;
		}
		
		par_true    = beta\mu\betap\mup\betan\mun;
		if (mod == "CNOPC") {
			par_true    = par_true\ron\rop
		}


	}else if (dgp == "NOP" || dgp == "NOPC"){ /* NOP */

		if (overlap == 1){  //@ no overlap @
			beta    = (.6, .4, .8 )';
			mu      = (0.26, 2.14)';
			betan   = (.2, .8, .3 )';
			betap   = (.4, .3, .9 )';
			mun     = -0.54;
			mup     =  0.54;

			xv   = 1 \ 0 \ 0 ;   
			znv  = 0 \ 1 \ 0 ;
			zpv  = 0 \ 0 \ 1 ;

		} else if (overlap == 2){ //  @ partial overlap @

			beta    = (.6, .4, .8 )';
			mu      = (0.21, 2.19)';
			betan   = (.2, .8, .3 )';
			betap   = (.4, .3, .9 )';
			mun     = -0.17 ;
			mup     = 0.68;

			xv      = 1 \ 1 \ 0 ;   
			znv     = 1 \ 0 \ 1 ;
			zpv     = 0 \ 1 \ 1 ;

		} else if (overlap == 3){  //@ complete overlap @
			xv   = 1 \ 1 \ 1 ;   
			znv  = 1 \ 1 \ 1 ;
			zpv  = 1 \ 1 \ 1 ;

			beta    = (.6, .4, .8 )';
			mu      = (0.09, 2.32)';
			betan   = (.2, .8, .3 )';
			betap   = (.4, .3, .9 )';
			mun     = -.72;
			mup     =  2.12;

		}

		beta        = select(beta,  xv);
		betan       = select(betan, znv);
		betap       = select(betap, zpv);

		par_true    = beta\mu\betap\mup\betan\mun;
		
		if (dgp=="NOPC") {
			rop = 0.6;
			ron = -0.1;
		} else {
			rop = 0;
			ron = 0;
		}
		if (mod == "NOPC") {
			par_true    = par_true\ron\rop
		}

	} else if(dgp == "OP"){ /* OP */

		beta    = (.6, .4, .8)';
		mu      = (-1.83, -1.01,  1.01,  1.83)';

		xv     = 0 \ 1 \ 1  ; 
		znv    = 0 \ 0 \ 0  ;
		zpv    = 0 \ 0 \ 0  ;
	  
		beta       = select(beta, xv);
		par_true   = beta\mu;
	}


	// names for printing (as for CNOP)
	{
		datanames = ("x1", "x2", "x3")'
		mu1names = "mu_" :+ strofreal(1::(cJ-1))

		bnames = "mu_" :+ strofreal(1::2)

		xnames = select(datanames, xv)
		zpnames = select(datanames, zpv)
		znnames = select(datanames, znv)

		if (mod == "CNOP" || mod == "CNOPC") {
			mu3namesp = "mup_" :+ strofreal(1::(cJ-inf))
			mu3namesn = "mun_" :+ strofreal(1::(inf-1))
			label  = xnames \ bnames \ zpnames \ mu3namesp \ znnames \ mu3namesn
		}
		if (mod == "NOP" || mod == "NOPC") {
			mu3namesp = "mup_" :+ strofreal(1::(cJ-inf-1))
			mu3namesn = "mun_" :+ strofreal(1::(inf-2))
			label  = xnames \ bnames \ zpnames \ mu3namesp \ znnames \ mu3namesn
		}
		if (mod == "OP") {
			label = label \ "rho1" \ "rho2"
		}
		
		if (mod == "CNOPC" || mod == "NOPC") {
			label = label \ "rho1" \ "rho2"
		}
	}
	/*******************************   END of MANUAL INPUT   ********************************/
	// lines 200-223 
	// some strange manipulations with data labels
	index = (xv+zpv+znv :> 0)'
	kxz = 3
	kx = sum(xv)
	kzp = sum(zpv)
	kzn = sum(znv)

	allvars = (x1,x2,x3)

	xdata = select(allvars, xv')
	zpdata = select(allvars, zpv')
	zndata = select(allvars, znv')

	// calculate values to compute ME at
	if (xbar==1) {
		xbar = select((2,0,0), index)
	} else if(xbar==0) {
		xbar = J(1,kxz,0)
	} else {
		xbar = select(xbar, index)
	}
	// calculate true ME
	me_true = asarray_create("real", 1, meNo)
	me0_true = asarray_create("real", 1, meNo)
	
	p_true = asarray_create("real", 1, meNo)
	p0_true = asarray_create("real", 1, meNo)
	
	class CNOPModel scalar true_model 
	true_model.model_class 	= mod
	true_model.params 	= par_true
	
	true_model.ncat 	= cJ
	true_model.infcat 	= inf
	true_model.corresp	= (xv' \ zpv' \ znv')
	
	for(i = 1; i<= meNo; i++) {
		generalME(par_true', allvars[i,], true_model,  (0,0,1), 1, _returnedME = .)
		asarray(me_true,i, rowshape(_returnedME,3))
		
		asarray(p_true,i, generalPredict(allvars[i,], true_model, 1))
		
		// INCOMOLETE: check LOOP argument
		if (mod != "OP" && mod != "NOP" && mod != "NOPC") {
			generalME(par_true', allvars[i,], true_model,  (0,0,1), 2, _returnedME = .)
			asarray(me0_true,i, rowshape(_returnedME,3))
			asarray(p0_true,i, generalPredict(allvars[i,], true_model, 2))
		}
	}

	// initialize counters
	{
	failure      = J(rep,1,0); // indicator of failures
	ret          = J(rep,1,0); // types of failures
	est_time     = J(rep,1,0); // in MS
	LL           = J(rep,1,0); // log likelihood

	b_op         = J(rep,kxz+cJ-1,0);
	bse_op       = J(rep,kxz+cJ-1,0);

	ME1_op       = J(rep,kxz*cJ,0); // had to fit to 2 dimensions!
	ME1_opse     = J(rep,kxz*cJ,0);
	CPme1_op     = J(rep,kxz*cJ,0);
	PRMSE_op     = J(rep,cJ,0);
	

	if (dgp == "OP") {

		b_nop        = J(rep,kxz+2+kxz+inf-2+kxz+cJ-inf-1,0);
		bse_nop      = J(rep,kxz+2+kxz+inf-2+kxz+cJ-inf-1,0);
		/*V_nop        = arrayinit(rep|kxz+2+kxz+inf-2+kxz+cJ-inf-1|kxz+2+kxz+inf-2+kxz+cJ-inf-1,0);*/
		b_tiop       = J(rep,kxz+2+kxz+inf-1+kxz+cJ-inf,0);
		bse_tiop     = J(rep,kxz+2+kxz+inf-1+kxz+cJ-inf,0);
		/*V_tiop       = arrayinit(rep|kxz+2+kxz+inf-1+kxz+cJ-inf|kxz+2+kxz+inf-1+kxz+cJ-inf,0);*/
		nargNOP		 = kxz+2+kxz+inf-2+kxz+cJ-inf-1
		nargCNOP	 = kxz+2+kxz+inf-1+kxz+cJ-inf
	} else {

		b_nop        = J(rep,kx+2+kzn+inf-2+kzp+cJ-inf-1,0);
		bse_nop      = J(rep,kx+2+kzn+inf-2+kzp+cJ-inf-1,0);
		/*V_nop        = arrayinit(rep|kx+2+kzn+inf-2+kzp+cJ-inf-1|kx+2+kzn+inf-2+kzp+cJ-inf-1,0);*/
		b_tiop       = J(rep,kx+2+kzn+inf-1+kzp+cJ-inf,0);
		bse_tiop     = J(rep,kx+2+kzn+inf-1+kzp+cJ-inf,0);
		nargNOP		 = kx+2+kzn+inf-2+kzp+cJ-inf-1
		nargCNOP 	 = kx+2+kzn+inf-1+kzp+cJ-inf
		/*V_tiop       = arrayinit(rep|kx+2+kzn+inf-1+kzp+cJ-inf|kx+2+kzn+inf-1+kzp+cJ-inf,0);*/
	}
	
	if (mod == "CNOP") {
		narg = nargCNOP
	} else if (mod == "CNOPC") {
		narg = nargCNOP + 2
	} else if (mod == "NOP") {
		narg = nargNOP
	} else if (mod == "NOPC") {
		narg = nargNOP + 2
	} else if (mod == "MIOPR") {
		narg = -1 // to think whether I need miopr
	} 

	
	b = J(rep, narg,0)
	bse = J(rep, narg,0)
	
	
	
	ME1_nop      = J(rep,kxz*cJ,0);
	ME1_nopse    = J(rep,kxz*cJ,0);
	CPme1_nop    = J(rep,kxz*cJ,0);
	PRMSE_nop    = J(rep,cJ,0);

	ME1_tiop     = J(rep,kxz*cJ,0);
	ME1_tiopse   = J(rep,kxz*cJ,0);
	CPme1_tiop   = J(rep,kxz*cJ,0);
	PRMSE_tiop   = J(rep,cJ,0);

	if (meNo > 0){
		all_me = J(rep*meNo,kxz*cJ,0);
		all_mese = J(rep*meNo,kxz*cJ,0);

		me_op        = J(rep*meNo,kxz*cJ,0);
		me_opse      = J(rep*meNo,kxz*cJ,0);
		/*CPme_op      = arrayinit(rep|meNo|kxz|cJ,0);*/
		me_nop       = J(rep*meNo,kxz*cJ,0);
		me_nopse     = J(rep*meNo,kxz*cJ,0);
		/*CPme_nop     = arrayinit(rep|meNo|kxz|cJ,0);*/
		me_tiop      = J(rep*meNo,kxz*cJ,0);
		me_tiopse    = J(rep*meNo,kxz*cJ,0);
		/*CPme_tiop    = arrayinit(rep|meNo|kxz|cJ,0);*/

	} else {
		all_me = J(rep,kxz*cJ,0);
		all_mese = J(rep,kxz*cJ,0);

		me_op        = J(rep,kxz*cJ,0);
		me_opse      = J(rep,kxz*cJ,0);
		/*CPme_op      = arrayinit(rep|1|kxz|cJ,0);*/
		me_nop       = J(rep,kxz*cJ,0);
		me_nopse     = J(rep,kxz*cJ,0);
		/*CPme_nop     = arrayinit(rep|1|kxz|cJ,0);*/
		me_tiop      = J(rep,kxz*cJ,0);
		me_tiopse    = J(rep,kxz*cJ,0);
		/*CPme_tiop    = arrayinit(rep|1|kxz|cJ,0);*/
		nmeCNOP		 = kxz*cJ
	}
	
	//me = J(rep, nmeCNOP,0)
	//mese = J(rep, nmeCNOP,0)
	
	probse = asarray_create("real", 2, rep)
	probse0 = asarray_create("real", 2, rep)
	
	mese = asarray_create("real", 2, rep)
	// 2d array, each element is 2d matrix, upper half is me, lower half is se
	mese0 = asarray_create("real", 2, rep)
	
	}
	
	//loops
	
	converged  	= J(rep, 1, 0)
	error_codes	= J(rep, 1, 0)
	
	for (jrep = 1; jrep <= rep; jrep++) {
		starttime = clock(c("current_time"),"hms")
		// switch(dgp)  create epsilons and y
		
		"ITERATION NUMBER: " + strofreal(jrep)

		
		if (dgp == "CNOP" || dgp == "CNOPC") { /* TIOP */
		
			eps  = rnormal(No,1,0,1);
			epsp = eps * rop + rnormal (No,1,0,1) * sqrt(1 - rop^2);
			epsn = eps * ron + rnormal (No,1,0,1) * sqrt(1 - ron^2);
			
			rs      = xdata*beta + eps;
			r       = J(No,3,0);
			
			r[,1]  = rs :< mu[1];
			r[,3]  = rs :>= mu[2];
			r[,2]  = 1 :- r[,1] :- r[,3];
			
			ysn     = zndata*betan + epsn;
			ysp     = zpdata*betap + epsp;
			
			yn                = J(No,rows(mun)+1,0);
			yn[,1]           = ysn :< mun[1];
			yn[,rows(mun)+1] = ysn :>= mun[rows(mun)];
			
			if (rows(mun) > 1){
				for (i = 2; i <=rows(mun); i++) {
					yn[.,i]   = (ysn :>= mun[i-1]) :* (ysn :< mun[i]);
				}
			}
			
			yp                = J(No,rows(mup)+1,0);
			yp[,1]           = ysp :< mup[1];
			yp[,rows(mup)+1] = ysp :>= mup[rows(mup)];
			
			if (rows(mup) > 1){
				for (i = 2; i <=rows(mup); i++) {
					yp[.,i]   = (ysp :>= mup[i-1]) :* (ysp :< mup[i]);
				}
			}
			
			yn      = yn* (-rows(mun)::0);
			yp      = yp*(0::(rows(mup)));
			
			yr      = yn, J(No,1,0),yp;
			y       = rowsum(r :* yr); /* observed y */
		
		} else if (dgp == "NOP" || dgp == "NOPC") {; /* NOP */
		
			eps  = rnormal(No,1,0,1);
			epsp = eps * rop + rnormal (No,1,0,1) * sqrt(1 - rop^2);
			epsn = eps * ron + rnormal (No,1,0,1) * sqrt(1 - ron^2);
			
			rs      = xdata*beta + eps;
			r       = J(No,3,0);
			
			r[,1]  = rs :< mu[1];
			r[,3]  = rs :>= mu[2];
			r[,2]  = 1 :- r[,1] :- r[,3];
			
			ysn     = zndata*betan + epsn;
			ysp     = zpdata*betap + epsp;
			
			yn                = J(No,rows(mun)+1,0);
			yn[,1]           = ysn :< mun[1];
			yn[,rows(mun)+1] = ysn :>= mun[rows(mun)];
			
			if (rows(mun) > 1){
				for (i =2; i<= rows(mun); i++){
					yn[,i]   = (ysn :>= mun[i-1]) :* (ysn :< mun[i]);
				}
			}
			
			yp                = J(No,rows(mup)+1,0);
			yp[,1]           = ysp :< mup[1];
			yp[,rows(mup)+1] = ysp :>= mup[rows(mup)];
			
			if (rows(mup) > 1){
				for (i =2; i<= rows(mup); i++){
					yp[,i]   = (ysp :>= mup[i-1]) :* (ysp :< mup[i]);
				}
			}
			
			yn      = yn*((-rows(mun)-1)::-1);
			yp      = yp*(1::(rows(mup)+1));
			
			yr      = yn,J(No,1,0),yp;
			y       = rowsum(r :* yr); /* observed y */
		
		} else if (dgp == "OP") { /* OP */
			
			eps  = rnormal(No,1,0,1);
			
			ys      = xdata*beta + eps;
			
			y       = J(No,cJ,0);
			y[,1]  = ys :< mu[1];
			y[,cJ] = ys :>= mu[cJ-1];
			
			if (cJ > 2){
				for (i = 2 ; i<=cJ-1; i++){
					y[,i]   = (ys :>= mu[i-1]) :* (ys :< mu[i]);
				}
			}
			
			y       = y*((1-inf)::(cJ-inf)); /* observed y */
		} else {
			"No such dgp implemented: " + dgp 
			y = .
		}
		// ESTIMATE EVERYTHING
		
		// WHAT are regressors exactly? should we use all data?
		// It's a bit unfair to allow OP to use all variables on the 1'st stage, and to restrict NOP and CNOP
		
		class CNOPModel scalar model
		
		if (mod == "CNOP") {
			model = estimateCNOP(y, xdata, zpdata, zndata, 0, quiet>0, ., robust, who)
		} else if (mod == "CNOPC") {
			model = estimateCNOPC(y, xdata, zpdata, zndata, 0, quiet>0, ., robust, who)
		} else if (mod == "NOP") {
			model = estimateNOP(y, xdata, zpdata, zndata, 0, quiet>0, ., robust, who)
		} else if (mod == "NOPC") {
			model = estimateNOPC(y, xdata, zpdata, zndata, 0, quiet>0, ., robust, who)
		} else if (mod == "MIOPR") {
			model = estimateMIOPR(y, xdata, zpdata, zndata, 0, quiet>0, ., robust, who)
		} 
		"01"
		if (model.error_code != 0) {
			est_time[jrep] = (clock(c("current_time"),"hms") - starttime)/1000
			converged[jrep] = 0
			error_codes[jrep] = model.error_code
			continue
		}
		"02"
		if (quiet < 2) {
			("", "coef", "se") \ (label, strofreal(model.params), strofreal(model.se))
		}
		"03"
		// flags that variables are the same
		corresp = (xv' \ zpv' \ znv')
		
		setcorresp(model, corresp)

		meandse = generalMEwithSE(xbar, model, dummiesVector, 1)
		"04"
		for(i = 1; i <= meNo; i++) {
			meandse = generalMEwithSE(allvars[i,], model, dummiesVector, 1)
			asarray(mese, (jrep, i), meandse)
			//meandse0 = generalMEwithSE(allvars[i,], model, dummiesVector, 2)
			//asarray(mese, (jrep, i), meandse0)
		}
		"05"
		// INCOMPLETE: remember estimation time and successfullness
		// INCOMPLETE: may need ME0
	
		
		
		/*
		yvalues = uniqrows(y)'
		prob1 = generalPredict(allvars,model1, 1)
		pred1 = rowsum((prob1:==rowmax(prob1)) :* yvalues)
		
		prob2 = generalPredict(allvars,model2, 1)
		pred2 = rowsum((prob2:==rowmax(prob2)) :* yvalues)
		*/
		"06"
		prob = generalPredict(allvars, model, 1)
		"07"
		for(i = 1; i <= meNo; i++) {
			pandse = generalPredictWithSE(allvars[i,], model, 1)
			asarray(probse, (jrep, i), pandse)
			//pandse0 = generalPredictWithSE(allvars[i,], model, 2)
			//asarray(probse0, (jrep, i), pandse0)
		}
		"08"
		// prob0 = generalPredict(allvars,model, 2)
		//pred = rowsum((prob3:==rowmax(prob)) :* yvalues)
		
		//rmse1 = sqrt(mean((pred1 - y):^2))
		//rmse2 = sqrt(mean((pred2 - y):^2))
		//rmse = sqrt(mean((pred3 - y):^2))
		
		//hit1 = mean(pred1:==y)
		//hit2 = mean(pred2:==y)
		//hit3 = mean(pred:==y)
		/*
		if (meNo > 0){
			mese = generalMEwithSE(
			all_me = J(rep*meNo,kxz*cJ,0);
			all_mese = J(rep*meNo,kxz*cJ,0)

		} else {
			all_me = J(rep,kxz*cJ,0);
			all_mese = J(rep,kxz*cJ,0);
		}
		*/
		
		// remember everything
		b[jrep, ] = model.params'
		bse[jrep, ] = model.se'
		//Hit[jrep,1] = hit3
		LL[rep] = model.logLik
		
		est_time[jrep] = (clock(c("current_time"), "hms") - starttime)/1000
		converged[jrep] = model.converged
		error_codes[jrep] = model.error_code
	} 
	
	//

	// save details

	// rearrange data to fit "ret" // whether iteration succeeded

	/***************************** DISPLAY RESULTS *********************************************/
	
	// which models to output
	filter = converged
	if (converged_only == 0) {
		filter = converged * 0 + 1
	}
	precise = converged
	if (precise_only > 0) {
		for (jrep = 1; jrep <= rep; jrep++) {
			if ((max(abs(b[jrep,]-par_true')) > precise_only) || max(bse[jrep,])> precise_only) {
				precise[jrep] = 0
				filter[jrep] = 0
			}
		}
	}
	frep = sum(filter)
	
	
	
	"ESTIMATION RESULTS"
	"Observations: " + strofreal(No) + ", estimations: " + strofreal(rep)
	"DGP: " + dgp + ", estimated model: " + mod +  ", overlap: " + strofreal(overlap)
	
	
	"Time: " + strofreal(sum(est_time)) + " seconds, on average " + strofreal(mean(est_time)) + " sec. per model, max " + strofreal(max(est_time)) + " sec. per model"
	"Share of converged models: " + strofreal(mean(converged))
	"Share of precise models: " + strofreal(mean(precise))
	"Share of critical errors: " + strofreal(mean(error_codes :!= 0))
	
	
	"Mean parameters:"

	meanparams = colsum(b:*filter)' :/ frep
	rmse =(rowsum(select((b' :- par_true):^2, filter')) :/ frep):^(0.5)
	
	meanse = colsum(bse:*filter)' :/ frep
	medianse = colMedians(select(bse, filter))'
	realse = (rowsum((b' :-meanparams):^2) :/ frep):^(0.5)

	cil = b - cv:*bse
	ciu = b + cv:*bse

	coverage = (par_true :> cil') :* (par_true :< ciu')
	meancoverage = rowsum(coverage:*filter') :/ frep

	compresult = ("", "true", "mean", "mean std", "real std", "mean/real std", "median/real std", "rmse", "coverage") \ (label, strofreal(par_true), strofreal(meanparams), strofreal(meanse), strofreal(realse),  strofreal(realse :/ meanse), strofreal(realse :/ medianse), strofreal(rmse), strofreal(meancoverage))
	compresult
	
	/********************** MARGINAL EFFECTS ***********************************/
	
	if (me_detailed_idx > 0 && me_detailed_idx <= meNo) {
		"Marginal effects: details at the "+strofreal(me_detailed_idx)+"th observation"
		trueme1 = asarray(me_true, me_detailed_idx)
		meanme1 = trueme1 * 0
		meanmese1 = trueme1 * 0
		meanmesq1 = trueme1 * 0
		mecp1 = trueme1 * 0
		mermse1 = trueme1 * 0
		memape1 = trueme1 * 0


		for(jrep = 1; jrep <=rep; jrep ++) {
			if (filter[jrep]) {
				meandse = asarray(mese, (jrep, me_detailed_idx))
				nrow = rows(meandse)
				me = meandse[1::(nrow/2), ]
				mse = meandse[((nrow/2)+1)::nrow, ]
				
				meanme1 = meanme1 + me/frep
				meanmese1 = meanmese1 + mse/frep
				meanmesq1 = meanmesq1 + (me:^2)/frep
				mecp1 = mecp1 + (trueme1 :> (me - cv*mse)) :* (trueme1 :< (me + cv*mse)) / frep
				mermse1 = mermse1 + ((trueme1 - me):^2)/frep
				memape1 = memape1 + abs(trueme1 - me)/frep
			}
		}
		truemese1 = sqrt(meanmesq1 - meanme1:^2)

		mermse1 = sqrt(mermse1)
		"true"
		trueme1
		"mean"
		meanme1
		"mean std"
		meanmese1
		"real std"
		truemese1
		"mean/real std"
		meanmese1 :/ truemese1
		"rmse"
		mermse1
		"mape"
		memape1
		"coverage"
		mecp1
	}


	
	if (meNo > 0 && meNo <= No) {
		"Marginal effects: aggregate results at the observations from 1st to "+strofreal(meNo)+"th"
		trueme3 = asarray(me_true, 1)
		meanme3 = trueme3 * 0
		meanmese3 = trueme3 * 0
		meanmesq3 = trueme3 * 0
		mecp3 = trueme3 * 0
		mermse3= trueme3 * 0
		memape3 = trueme3 * 0
		truemese3= trueme3 * 0
		meseratio3 = trueme3 * 0


		for(i = 1; i <= meNo; i++) {
			trueme2 = asarray(me_true, i)
			meanme2 = trueme2 * 0
			meanmese2 = trueme2 * 0
			meanmesq2 = trueme2 * 0
			mecp2 = trueme2 * 0
			mermse2 = trueme2 * 0
			memape2 = trueme2 * 0
			for(jrep = 1; jrep <=rep; jrep ++) {
				if (filter[jrep]) {
					meandse = asarray(mese, (jrep, i))
					nrow = rows(meandse)
					me = meandse[1::(nrow/2), ]
					mse = meandse[((nrow/2)+1)::nrow, ]
					
					meanme2 = meanme2 + me/frep
					meanmese2 = meanmese2 + mse/frep
					meanmesq2 = meanmesq2 + (me:^2)/frep
					mecp2 = mecp2 + (trueme2 :> (me - cv*mse)) :* (trueme2 :< (me + cv*mse)) / frep
					mermse2 = mermse2 + ((trueme2 - me):^2)/frep
					memape2 = memape2 + abs(trueme2 - me)/frep
				}
			}
			truemese2 = sqrt(meanmesq2 - meanme2:^2)
			meseratio2 = (meanmese2 :/ truemese2)
			
			
			trueme3 	= trueme3 + trueme2/meNo
			meanme3 	= meanme3 + meanme2/meNo
			meanmese3 	= meanmese3 + meanmese2/meNo
			truemese3 	= truemese3 + truemese2/meNo
			mermse3		= mermse3 + mermse2/meNo
			meseratio3 	= meseratio3 + meseratio2/meNo
			memape3		= memape3 + memape2/meNo
			mecp3		= mecp3 + mecp2/meNo
		}

		"Average RMSE of Marginal Effects"
		mean(mean(mermse3)')
		"Average MAPE of Marginal Effects"
		mean(mean(memape3)')
		"Average Coverage Probabilities of MEs:"
		mean(mean(mecp3)')
		"Average ratio of estimated st. errors of MEs and st. deviation of estimated MEs"
		mean(mean(meseratio3)')
	}
	
	/********************** PROBABILITIES ***********************************/
	if (meNo > 0 && meNo <= No) {
		"Probabilities: aggregate results at the observations from 1st to "+strofreal(meNo)+"th"
		truep3 = asarray(p_true, 1)
		meanp3 = truep3 * 0
		meanpse3 = truep3 * 0
		meanpsq3 = truep3 * 0
		pcp3 = truep3 * 0
		prmse3= truep3 * 0
		pmape3 = truep3 * 0
		truepse3= truep3 * 0
		pseratio3 = truep3 * 0


		for(i = 1; i <= meNo; i++) {
			truep2 = asarray(p_true, i)
			meanp2 = truep2 * 0
			meanpse2 = truep2 * 0
			meanpsq2 = truep2 * 0
			pcp2 = truep2 * 0
			prmse2 = truep2 * 0
			pmape2 = truep2 * 0
			for(jrep = 1; jrep <=rep; jrep ++) {
				if (filter[jrep]) {
					pandse = asarray(probse, (jrep, i))
					p = pandse[1, ]
					pse = pandse[2, ]
					
					meanp2 = meanp2 + p/frep
					meanpse2 = meanpse2 + pse/frep
					meanpsq2 = meanpsq2 + (p:^2)/frep
					pcp2 = pcp2 + (truep2 :> (p- cv*pse)) :* (truep2 :< (p + cv*pse)) / frep
					prmse2 = prmse2 + ((truep2 - p):^2)/frep
					pmape2 = pmape2 + abs(truep2 - p)/frep
				}
			}
			truepse2 = sqrt(meanpsq2 - meanp2:^2)
			pseratio2 = (meanpse2 :/ truepse2)
			
			
			truep3 	= truep3 + truep2/meNo
			meanmp3 	= meanp3 + meanp2/meNo
			meanpse3 	= meanpse3 + meanpse2/meNo
			truepse3 	= truepse3 + truepse2/meNo
			prmse3		= prmse3 + prmse2/meNo
			pseratio3 	= pseratio3 + pseratio2/meNo
			pmape3		= pmape3 + pmape2/meNo
			pcp3		= pcp3 + pcp2/meNo
		}

		"Average RMSE of Probabilities"
		mean(mean(prmse3)')
		"Average MAPE of Probabilities"
		mean(mean(pmape3)')
		"Average Coverage Probabilities of Probabilities:"
		mean(mean(pcp3)')
		"Average ratio of estimated st. errors of Probabilities and st. deviation of estimated Probabilities"
		mean(mean(pseratio3)')
	}
	// 
	if (filesuffix != .) {
		unlink("par_" + filesuffix)
		unlink("mes_sum_" + filesuffix)
		unlink("prob_sum_" + filesuffix)
		fh  = fopen("par_" + filesuffix, "w")
		fputmatrix(fh,(dgp, mod, strofreal(No), strofreal(rep), strofreal(meNo), strofreal(overlap), strofreal(sum(est_time)), strofreal(mean(converged))))
		fputmatrix(fh,converged)
		fputmatrix(fh,par_true)
		fputmatrix(fh,compresult)
		fputmatrix(fh,b)
		fputmatrix(fh,bse)
		fclose(fh)
		if (meNo > 0 && meNo <= No) {
			fh  = fopen("mes_sum_" + filesuffix, "w")
			fputmatrix(fh, mermse3)
			fputmatrix(fh, meseratio3)
			fputmatrix(fh, memape3)
			fputmatrix(fh, mecp3)
			fclose(fh)
			fh  = fopen("prob_sum_" + filesuffix, "w")
			fputmatrix(fh, prmse3)
			fputmatrix(fh, pseratio3)
			fputmatrix(fh, pmape3)
			fputmatrix(fh, pcp3)
			fclose(fh)	
		}
	}
	//
	
}

// And here I run it

No = 200
rep = 5000 
meNo = 1
me_detailed_idx = 0 // value between 1 and meNo if want details
converged_only = 1  // exclude models which did not converge to maximum of likelihood
precise_only = 5  // error of b or value of bse, starting with which model is considered as bad convergence

// all with 
nos = (200, 500, 1000)
mods = ("CNOPC")
overlaps = (1,2,3)

rseed(42)
doSimulations(filesuffix = "tmp_", No, rep, meNo, mod = "CNOP", dgp = "CNOP", overlap = 2, allvars = x[1..1000,1..3], cv = invnormal(.975), quiet = 2, seed = 1, me_detailed_idx, converged_only, precise_only, robust = 0, who = .,  b = ., bse=., mese=., mese0 =., probse=., probse0 = ., par_true = ., me_true = ., me0_true = ., p_true = ., p0_true = ., label = ., est_time = ., converged = ., error_codes = ., precise = .)


/*
for(k=1; k <=4; k++) {
	for(j = 1; j <= 3; j++) {
		for (i = 1; i <= 3; i++) {		
			hash = strofreal(nos[i]) + "_" + strofreal(overlaps[j]) + "_" + mods[k]
			doSimulations(filesuffix = "tst_" + hash, nos[i], rep, meNo, mod = mods[k], dgp = mods[k], overlap = overlaps[j], allvars = x[1..1000,1..3], cv = invnormal(.975), quiet = 2, seed = 1, me_detailed_idx, converged_only, robust = 0, who = .,  b = ., bse=., mese=., mese0 =., probse=., probse0 = ., par_true = ., me_true = ., me0_true = ., p_true = ., p0_true = ., label = ., est_time = ., converged = ., error_codes = .)
		}
	}
}
*/
end










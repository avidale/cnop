version 12

mata

class CNOPModel {
	string scalar model_class
	
	// utility parameters
	real vector beta
	real vector a
	real vector gamma
	real vector mu
	real scalar n
	real scalar k
	
	real scalar ncat
	real scalar infcat
	real vector allcat
	real scalar ncatp
	real scalar ncatn
	
	// parameters for positive-negative
	real scalar kp
	real scalar kn
	real vector gammap
	real vector mup
	real vector gamman
	real vector mun
	real scalar rop
	real scalar ron
	
	// estimation options
	real scalar robust
	
	// main estimation results
	real vector params
	real matrix V
	real matrix V_rob
	real vector se
	real vector t
	real vector se_rob
	real vector t_rob
	real scalar logLik
	real vector pval
	
	// probabilities
	real matrix probabilities
	// probability of r == 1
	real scalar p1_sorted
	real scalar p1_range
	real scalar p1_xbar
	
	// optimization outcome
	string scalar retCode 
	real scalar etime
	real scalar error_code
	real scalar converged
	real scalar iterations
	
	// infocriteria
	real scalar AIC		
	real scalar BIC		
	real scalar CAIC	
	real scalar AICc	
	real scalar HQIC
	real scalar R2 // mcFadden R2 = 1 - lnL(model)/lnL(simple model)
	real scalar logLik0
	
	// marginal effects
	real vector me
	real vector mese
	real vector met
	real vector mepval
	
	// description of data
	string scalar yname
	string scalar xnames
	string scalar znames
	string scalar zpnames
	string scalar znnames
	
	string vector XZnames	// Stata names of independent variables
	real matrix corresp	// coincidence of independent variables in X, Zp and Zn matrices
	/*
	Example of corresp:
	
	*/
	real vector XZmeans		// means of independent variables
	string vector eqnames
	string vector parnames
}


function Vuongtest(class CNOPModel scalar model1, class CNOPModel scalar model2, real matrix trueClassification){
	if(model1.n != model2.n){
		 _error("Vuong test failed: two models are estimated on different numpber of observations")
	}
	if(model1.ncat != model2.ncat){
		 _error("Vuong test failed: two models include different number of categories")
	}
	f1      = colsum((model1.probabilities :* trueClassification)');
	f2      = colsum((model2.probabilities :* trueClassification)');
	mi      = ln(f1 :/ f2);
	vtop    = mean(mi') * sqrt(model1.n);
	vbot    = sqrt(variance(mi'));
	return(vtop / vbot)
}

function Vuong(real matrix prob1, real matrix prob2, real matrix trueClassification){
	f1      = colsum((prob1 :* trueClassification)');
	f2      = colsum((prob2 :* trueClassification)');
	mi      = ln(f1 :/ f2);
	vtop    = mean(mi') * sqrt(rows(prob1));
	vbot    = sqrt(variance(mi'));
	return(vtop / vbot)
}
  
end

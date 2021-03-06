version 12

clear matrix

mata: mata clear
// subroutines from .ado files will be saved in MATA libraries upon completion

//base: C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\

// define CNOPModel structure
run CNOPishModel_definition.ado
// define all auxiliary functions for optimization and ME
run gradients.ado
// define functions for model estimation
run inflatedOP_estimation_routines.ado

capture program drop ziop3 
capture program drop nop 
capture program drop ziopmargins
capture program drop ziopprobabilities
capture program drop ziopcontrasts
capture program drop ziop2
capture program drop zioppredict
capture program drop ziopclassification
capture program drop ziopvuong


// confusion matrix (classification table) for the last ziop-like command
program ziopclassification, rclass
	version 13
	syntax [if] [in]
	marksample touse
	// this is the generation part, specific to the ZIOP-like commands
	predict _predicted, output(mode)
	label variable _predicted "Predicted outcomes"
	gen _actual = `e(depvar)'
	label variable _actual "Actual outcomes"
	// this is the general comparison part
	gen _correct_predicted = _predicted == _actual
	display "Classification table"
	tab _actual _predicted if `touse', matcell(cells) matrow(labels)
	quietly sum _correct_predicted if `touse'
	mata: printf("Accuracy (%% of correct predictions) = %9.4f \n", `r(mean)') 
	mata: printf("Brier score                         = %9.4f \n", CNOP_last_model.brier_score)
	mata: printf("Ranked probability score            = %9.4f \n", CNOP_last_model.ranked_probability_score)
	display ""
	mata: classification_calc_large("_actual", "_predicted", "`touse'")
	drop _predicted _correct_predicted _actual
	return local accuracy = `r(mean)'
	return matrix noise result
end

// vuong test to compare two models
program ziopvuong, rclass
	version 13
	args modelspec1 modelspec2
	// todo: try to est store current environment
	quietly est restore `modelspec1'
	mat ll_1 = e(ll_obs)
	local k_1 = e(k)
	quietly est restore `modelspec2'
	mat ll_2 = e(ll_obs)
	local k_2 = e(k)
	mat ll_diff = ll_1 - ll_2
	display "Vuong non-nested test for `modelspec1' vs `modelspec2'"
	mata: vuong_calc()
	// todo: try to restore current environment
	return scalar mean_diff = mean_diff
	return scalar std_diff = std_diff
	return scalar N = n_obs
	return scalar vuong = vuong
	return scalar vuong_aic = vuongAIC
	return scalar vuong_bic = vuongBIC
	return scalar pvalue = pvalue
	return scalar pvalue_aic = pvalueAIC
	return scalar pvalue_bic = pvalueBIC
end


// prediction for NOP, ZIOP2 and ZIOP3
program zioppredict
	version 13
	syntax name [if] [in] [, zeros regimes output(string asis)]
	marksample touse
	mata: CNOP_predict(CNOP_last_model, "`1'", "`zeros'" == "zeros", "`regimes'"=="regimes", "`touse'", "`output'")
end

// produces marginal effects for NOP, ZIOP2 and ZIOP3
program ziopmargins, rclass
	version 13
	syntax [, at(string asis) zeros regimes]
	mata: CNOPmargins(CNOP_last_model, "`at'", "`zeros'" == "zeros", "`regimes'"=="regimes")
	return matrix at at
	return matrix me me
	return matrix se se
	return matrix t t
	return matrix pval pval
end

program ziopprobabilities, rclass
	version 13
	syntax [, at(string asis) zeros regimes]
	mata: CNOPprobabilities(CNOP_last_model, "`at'", "`zeros'" == "zeros", "`regimes'"=="regimes")
	return matrix at at
	return matrix me me
	return matrix se se
	return matrix t t
	return matrix pval pval
end

program ziopcontrasts, rclass
	version 13
	syntax [, at(string asis) to(string asis) zeros regimes]
	mata: CNOPcontrasts(CNOP_last_model, "`at'", "`to'", "`zeros'" == "zeros", "`regimes'"=="regimes")
	return matrix between between
	return matrix me me
	return matrix se se
	return matrix t t
	return matrix pval pval
end


// estimates ZIOP3 and ZIOP3(c) models
program ziop3, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, POSindepvars(varlist) NEGindepvars(varlist) INFcat(real 0) ENDOswitch CLUster(varname) ROBust INITial(string asis) nolog VUong]
	marksample touse
	mata: CNOP_last_model = processCNOP("`varlist'", "`posindepvars'", "`negindepvars'", `infcat', "`endoswitch'" == "endoswitch", "`touse'", "`robust'" == "robust", "`cluster'", "`initial'", "`log'" == "nolog")

	ereturn post b V, esample(`touse') obs(`N') depname(`depvar')
	ereturn local predict "zioppredict"
	ereturn local cmd "ziop3"
	ereturn scalar ll = ll
	ereturn scalar k = k
	ereturn matrix ll_obs ll_obs
	ereturn scalar r2_p = r2_p
	ereturn scalar k_cat = k_cat
	ereturn scalar df_m = df_m
	ereturn scalar ll_0 = ll_0
	ereturn scalar chi2 = chi2
	ereturn scalar p = p
	ereturn scalar aic = aic
	ereturn scalar bic = bic
	ereturn display
	if "`vuong'" == "vuong" {
		display "Vuong test versus ordered probit:"
		mata: vuong_vs_op(CNOP_last_model)
		ereturn scalar vuong = vuong
		ereturn scalar vuong_aic = vuongAIC
		ereturn scalar vuong_bic = vuongBIC
	}
end

// estimates MIOP(r) model
program ziop2, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, OUTindepvars(varlist) INFcat(real 0) ENDOswitch CLuster(varname) ROBust INITial(string asis) nolog]
	marksample touse
	mata: CNOP_last_model = processMIOPR("`varlist'", "`outindepvars'", `infcat', "`endoswitch'" == "endoswitch", "`touse'", "`robust'"=="robust","`cluster'", "`initial'", "`log'" == "nolog")
	
	ereturn post b V, esample(`touse')  depname(`depvar') obs(`N')
	ereturn local predict "zioppredict"
	ereturn local cmd "ziop2"
	ereturn scalar ll = ll
	ereturn scalar k = k
	ereturn matrix ll_obs ll_obs
	ereturn scalar r2_p = r2_p
	ereturn scalar k_cat = k_cat
	ereturn scalar df_m = df_m
	ereturn scalar ll_0 = ll_0
	ereturn scalar chi2 = chi2
	ereturn scalar p = p
	ereturn scalar aic = aic
	ereturn scalar bic = bic
	ereturn display
end

// estimates NOP and NOP(c) models
program nop, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, POSindepvars(varlist) NEGindepvars(varlist) INFcat(real 0) ENDOswitch CLuster(varname) ROBust INITial(string asis) nolog VUong]
	marksample touse
	mata: CNOP_last_model = processNOP("`varlist'", "`posindepvars'", "`negindepvars'", `infcat', "`endoswitch'" == "endoswitch", "`touse'", "`robust'" == "robust", "`cluster'", "`initial'", "`log'" == "nolog")

	ereturn post b V, esample(`touse') obs(`N') depname(`depvar')
	ereturn local predict "zioppredict"
	ereturn local cmd "nop"
	ereturn scalar ll = ll
	ereturn scalar k = k
	ereturn matrix ll_obs ll_obs
	ereturn scalar r2_p = r2_p
	ereturn scalar k_cat = k_cat
	ereturn scalar df_m = df_m
	ereturn scalar ll_0 = ll_0
	ereturn scalar chi2 = chi2
	ereturn scalar p = p
	ereturn scalar aic = aic
	ereturn scalar bic = bic
	ereturn display
	if "`vuong'" == "vuong" {
		display "Vuong test versus ordered probit:"
		mata: vuong_vs_op(CNOP_last_model)
		ereturn scalar vuong = vuong
		ereturn scalar vuong_aic = vuongAIC
		ereturn scalar vuong_bic = vuongBIC
	}
end


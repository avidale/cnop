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
capture program drop ZIOP_predict
capture program drop ziopconfusion


// confusion matrix (classification table) for the last ziop-like command
program ziopconfusion
	version 13
	predict _predicted, output(mode)
	gen _correct_predicted = _predicted == `e(depvar)'
	display "Classification table"
	tab _predicted `e(depvar)'
	quietly sum _correct_predicted
	display "% Correctly Predicted = " round(`r(mean)', 0.0001)
	drop _predicted _correct_predicted 
end

// prediction for NOP, ZIOP2 and ZIOP3
program ZIOP_predict
	version 13
	syntax name [if] [in] [, zeros regime output(string asis) at(string asis) point]
	// "at" and "point" are not implemented yet
	marksample touse
	mata: CNOP_predict(CNOP_last_model, "`1'", "`zeros'" == "zeros", "`regime'"=="regime", "`touse'", "`output'")
end

// produces marginal effects for NOP, ZIOP2 and ZIOP3
program ziopmargins
	version 13
	syntax [, at(string asis) nominal(varlist) zeros regime]
	mata: CNOPmargins(CNOP_last_model, "`at'", "`nominal'", "`zeros'" == "zeros", "`regime'"=="regime")
	display "Marginal effects of all variables on probabilities"
	mat list r(me)
	display "Standard errors of marginal effects"
	mat list r(se)
end

program ziopprobabilities
	version 13
	syntax [, at(string asis) zeros regime]
	mata: CNOPprobabilities(CNOP_last_model, "`at'", "`zeros'" == "zeros", "`regime'"=="regime")
	display "Predicted probabilities"
	mat list r(me)
	display "Standard errors of probabilities"
	mat list r(se)
end

program ziopcontrasts
	version 13
	syntax [, at(string asis) to(string asis) zeros regime]
	mata: CNOPcontrasts(CNOP_last_model, "`at'", "`to'", "`zeros'" == "zeros", "`regime'"=="regime")
	display "Contrasts of predicted probabilities"
	mat list r(me)
	display "Standard errors of contrasts"
	mat list r(se)
end


// estimates ZIOP3 and ZIOP3(c) models
program ziop3, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, xp(varlist) xn(varlist) infcat(integer 0) endoswitch cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processCNOP("`varlist'", "`xp'", "`xn'", `infcat', "`endoswitch'" == "endoswitch", "`touse'", "`robust'" == "robust", "`cluster'", "`initial'")

	ereturn post b V, esample(`touse') obs(`N') depname(`depvar')
	ereturn local predict "ZIOP_predict"
	ereturn local cmd "ziop3"
	ereturn display
end

// estimates MIOP(r) model
program ziop2, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, x(varlist) infcat(integer 0) endoswitch cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processMIOPR("`varlist'", "`x'", `infcat', "`endoswitch'" == "endoswitch", "`touse'", "`robust'"=="robust","`cluster'", "`initial'")
	
	ereturn post b V, esample(`touse')  depname(`depvar') obs(`N')
	ereturn local predict "ZIOP_predict"
	ereturn local cmd "ziop2"
	
	ereturn display
end

// estimates NOP and NOP(c) models
program nop, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, xp(varlist) xn(varlist) infcat(integer 0) endoswitch cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processNOP("`varlist'", "`xp'", "`xn'", `infcat', "`endoswitch'" == "endoswitch", "`touse'", "`robust'" == "robust", "`cluster'", "`initial'")

	ereturn post b V, esample(`touse') obs(`N') depname(`depvar')
	ereturn local predict "ZIOP_predict"
	ereturn local cmd "nop"
	ereturn display
end


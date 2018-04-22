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
capture program drop miop
capture program drop ZIOP_predict



// prediction for NOP, ZIOP2 and ZIOP3
program ZIOP_predict
	version 13
	syntax name [if] [in] [, zeroes regime output(string asis) at(string asis) point]
	// "at" and "point" are not implemented yet
	marksample touse
	mata: CNOP_predict(CNOP_last_model, "`1'", "`zeroes'" == "zeroes", "`regime'"=="regime", "`touse'", "`output'")
end

// produces marginal effects for NOP, ZIOP2 and ZIOP3
program ziopmargins
	version 13
	syntax [, at(string asis) nominal(varlist) zeroes regime]
	mata: CNOPmargins(CNOP_last_model, "`at'", "`nominal'", "`zeroes'" == "zeroes", "`regime'"=="regime")
	display "Marginal effects of all variables on probabilities"
	mat list r(me)
	display "Standard errors of marginal effects"
	mat list r(se)
end

program ziopprobabilities
	version 13
	syntax [, at(string asis) zeroes regime]
	mata: CNOPprobabilities(CNOP_last_model, "`at'", "`zeroes'" == "zeroes", "`regime'"=="regime")
	display "Predicted probabilities"
	mat list r(me)
	display "Standard errors of probabilities"
	mat list r(se)
end

program ziopcontrasts
	version 13
	syntax [, at(string asis) to(string asis) zeroes regime]
	mata: CNOPcontrasts(CNOP_last_model, "`at'", "`to'", "`zeroes'" == "zeroes", "`regime'"=="regime")
	display "Contrasts of predicted probabilities"
	mat list r(me)
	display "Standard errors of contrasts"
	mat list r(se)
end


// estimates ZIOP3 and ZIOP3(c) models
program ziop3, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, zp(varlist) zn(varlist) infcat(integer 0) correlated cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processCNOP("`varlist'", "`zp'", "`zn'", `infcat', "`correlated'" == "correlated", "`touse'", "`robust'" == "robust", "`cluster'", "`initial'")

	ereturn post b V, esample(`touse') obs(`N') depname(`depvar')
	ereturn local predict "ZIOP_predict"
	ereturn local cmd "ziop3"
	ereturn display
end

// estimates MIOP(r) model
program ziop2, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, z(varlist) infcat(integer 0) correlated cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processMIOPR("`varlist'", "`z'", `infcat', "`correlated'" == "correlated", "`touse'", "`robust'"=="robust","`cluster'", "`initial'")
	
	ereturn post b V, esample(`touse')  depname(`depvar') obs(`N')
	ereturn local predict "ZIOP_predict"
	ereturn local cmd "ziop2"
	
	ereturn display
end

// estimates NOP and NOP(c) models
program nop, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, zp(varlist) zn(varlist) infcat(integer 0) correlated cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processNOP("`varlist'", "`zp'", "`zn'", `infcat', "`correlated'" == "correlated", "`touse'", "`robust'" == "robust", "`cluster'", "`initial'")

	ereturn post b V, esample(`touse') obs(`N') depname(`depvar')
	ereturn local predict "ZIOP_predict"
	ereturn local cmd "nop"
	ereturn display
end


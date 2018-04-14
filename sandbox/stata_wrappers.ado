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

capture program drop cnop 
capture program drop nop 
capture program drop cnopmargins
capture program drop cnopprobabilities
capture program drop cnopcontrasts
capture program drop miop
capture program drop CNOP_predict



// prediction for MIOP(r), CNOP, CNOP(c)
program CNOP_predict
	version 13
	syntax name [if] [in] [, zeroes regime output(string asis) at(string asis) point]
	// "at" and "point" are not implemented yet
	marksample touse
	mata: CNOP_predict(CNOP_last_model, "`1'", "`zeroes'" == "zeroes", "`regime'"=="regime", "`touse'", "`output'")
end

// produces marginal effects for MIOP(r), CNOP, CNOP(c)
program cnopmargins
	version 13
	syntax [, at(string asis) nominal(varlist) zeroes regime]
	mata: CNOPmargins(CNOP_last_model, "`at'", "`nominal'", "`zeroes'" == "zeroes", "`regime'"=="regime")
	display "Marginal effects of all variables on probabilities"
	mat list r(me)
	display "Standard errors of marginal effects"
	mat list r(se)
end

program cnopprobabilities
	version 13
	syntax [, at(string asis) zeroes regime]
	mata: CNOPprobabilities(CNOP_last_model, "`at'", "`zeroes'" == "zeroes", "`regime'"=="regime")
	display "Predicted probabilities"
	mat list r(me)
	display "Standard errors of probabilities"
	mat list r(se)
end

program cnopcontrasts
	version 13
	syntax [, at(string asis) to(string asis) zeroes regime]
	mata: CNOPcontrasts(CNOP_last_model, "`at'", "`to'", "`zeroes'" == "zeroes", "`regime'"=="regime")
	display "Contrasts of predicted probabilities"
	mat list r(me)
	display "Standard errors of contrasts"
	mat list r(se)
end


// estimates CNOP and CNOP(c) models
program cnop, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, zp(varlist) zn(varlist) infcat(integer 0) correlated cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processCNOP("`varlist'", "`zp'", "`zn'", `infcat', "`correlated'" == "correlated", "`touse'", "`robust'" == "robust", "`cluster'", "`initial'")

	ereturn post b V, esample(`touse') obs(`N') depname(`depvar')
	ereturn local predict "CNOP_predict"
	ereturn local cmd "cnop"
	ereturn display
end

// estimates MIOP(r) model
program miop, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, z(varlist) infcat(integer 0) correlated cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processMIOPR("`varlist'", "`z'", `infcat', "`correlated'" == "correlated", "`touse'", "`robust'"=="robust","`cluster'", "`initial'")
	
	ereturn post b V, esample(`touse')  depname(`depvar') obs(`N')
	ereturn local predict "CNOP_predict"
	ereturn local cmd "miop"
	
	ereturn display
end

// estimates NOP and NOP(c) models
program nop, eclass
	version 13
	syntax varlist(min=2) [if] [in] [, zp(varlist) zn(varlist) infcat(integer 0) correlated cluster(varname) robust initial(string asis)]
	marksample touse
	mata: CNOP_last_model = processNOP("`varlist'", "`zp'", "`zn'", `infcat', "`correlated'" == "correlated", "`touse'", "`robust'" == "robust", "`cluster'", "`initial'")

	ereturn post b V, esample(`touse') obs(`N') depname(`depvar')
	ereturn local predict "CNOP_predict"
	ereturn local cmd "nop"
	ereturn display
end


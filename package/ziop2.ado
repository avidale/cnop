
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


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


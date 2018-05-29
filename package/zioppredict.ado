
// prediction for NOP, ZIOP2 and ZIOP3
program zioppredict
	version 13
	syntax name [if] [in] [, zeros regimes output(string asis)]
	marksample touse
	mata: CNOP_predict(CNOP_last_model, "`1'", "`zeros'" == "zeros", "`regimes'"=="regimes", "`touse'", "`output'")
end

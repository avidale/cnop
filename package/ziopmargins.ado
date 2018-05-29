
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

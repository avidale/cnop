
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

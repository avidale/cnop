
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

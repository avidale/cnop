
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

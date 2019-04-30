// confusion matrix (classification table) for the last ziop-like command
program ziopclassification, rclass
	version 13
	syntax [if] [in]
	marksample touse
	// this is the generation part, specific to the ZIOP-like commands
	predict _predicted, output(mode)
	label variable _predicted "Predicted outcomes"
	gen _actual = `e(depvar)'
	label variable _actual "Actual outcomes"
	// this is the general comparison part
	gen _correct_predicted = _predicted == _actual
	display "Classification table"
	tab _actual _predicted if `touse', matcell(cells) matrow(labels)
	quietly sum _correct_predicted if `touse'
	mata: printf("Accuracy (%% of correct predictions) = %9.4f \n", `r(mean)') 
	mata: printf("Brier score                         = %9.4f \n", CNOP_last_model.brier_score)
	mata: printf("Ranked probability score            = %9.4f \n", CNOP_last_model.ranked_probability_score)
	display ""
	mata: classification_calc_large("_actual", "_predicted", "`touse'")
	drop _predicted _correct_predicted _actual
	return local accuracy = `r(mean)'
	return matrix noise result
end

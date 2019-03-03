
mata

mid = colsum(all_crit1)' / rows(all_crit1), colsum(all_crit2)' / rows(all_crit2)
btm = colsum(all_crit12)' / rows(all_crit12)


fname = "MC results feb 2019.xlsx"
/*
sheetname = "MC results (10it_conv)"
*/

excel = xl()
excel.load_book(fname)
excel.set_mode("open")
excel.set_sheet(sheetname)

if(DGP == "OP") {
	left = 4
} else {
	left = 18
}
shift = 7
/*left = left + (repeat_dataset - 1) * 10*/

excel.put_number(8, left,       effp1_new)
excel.put_number(8, left+shift, effp2_new)
excel.put_number(12, left,       prmse1)
excel.put_number(12, left+shift, prmse2)
excel.put_number(15, left,       effm1_new)
excel.put_number(15, left+shift, effm2_new)


excel.put_number(32, left, colsum(all_cmp1)' / rows(all_cmp1))

excel.put_number(32, left+shift, colsum(all_cmp2)' / rows(all_cmp2))

if (MDLS[1] == DGP) {
	top = 49
} else if (MDLS[2] == DGP) {
	top = 55
} 

excel.put_number(top, 4, effx_new)

excel.put_number(61, left, (it, estimations, ready, con1, con2))

excel.put_number(64, left, mean_criteria)

excel.put_number(72, left, mean_comparison)

/*
excel.put_number(42, left, mid)
excel.put_number(50, left, btm)
*/

excel.close_book()

end

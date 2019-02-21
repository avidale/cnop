
mata

mid = colsum(all_crit1)' / rows(all_crit1), colsum(all_crit2)' / rows(all_crit2)
btm = colsum(all_crit12)' / rows(all_crit12)


fname = "MC results feb 2019.xlsx"
sheetname = "MC results (2)"

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

excel.put_number(7, left,       effp1_new)
excel.put_number(7, left+shift, effp2_new)
excel.put_number(14, left,       effm1_new)
excel.put_number(14, left+shift, effm2_new)


excel.put_number(31, left, colsum(all_cmp1)' / rows(all_cmp1))

excel.put_number(31, left+shift, colsum(all_cmp2)' / rows(all_cmp2))

/*
excel.put_number(42, left, mid)
excel.put_number(50, left, btm)
*/

excel.close_book()

end

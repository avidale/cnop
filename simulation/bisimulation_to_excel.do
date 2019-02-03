
mata

top = column1', column2'
mid = colsum(all_crit1)' / rows(all_crit1), colsum(all_crit2)' / rows(all_crit2)
btm = colsum(all_crit12)' / rows(all_crit12)


fname = "sim_compare_results.xlsx"
sheetname = "MC results (20190202)"

excel = xl()
excel.load_book(fname)
excel.set_mode("open")
excel.set_sheet(sheetname)

if(DGP == "OP") {
	left = 4
} else {
	left = 7
}

excel.put_number(7, left, top[(1,2,3,4),])
excel.put_number(12, left, top[(5,6,7,8),])
excel.put_number(17, left, top[(9,10,11,12,13,14,15,16,17,18,19,20,21),])

excel.put_number(34, left, mid)
excel.put_number(42, left, btm)

excel.close_book()

end

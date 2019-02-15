
mata

top = column1', column2'
mid = colsum(all_crit1)' / rows(all_crit1), colsum(all_crit2)' / rows(all_crit2)
btm = colsum(all_crit12)' / rows(all_crit12)


fname = "sim_compare_results.xlsx"
sheetname = "MC results (20190205)-10.0"

excel = xl()
excel.load_book(fname)
excel.set_mode("open")
excel.set_sheet(sheetname)

if(DGP == "OP") {
	left = 4
} else {
	left = 7
}
left = left + (repeat_dataset - 1) * 10

excel.put_number(7, left, top[range(1,12,1),])
excel.put_number(20, left, top[range(13,16,1),])
excel.put_number(25, left, top[range(17,29,1),])

excel.put_number(42, left, mid)
excel.put_number(50, left, btm)

excel.close_book()

end

program varsum
version 12
syntax varname [if] [in]
marksample touse
mata: calcsum("`varlist'", "`touse'")
display as txt " sum = " as res r(sum)
end
version 12
mata:
void calcsum(string scalar varname, string scalar touse)
{
real colvector x
st_view(x, ., varname, touse)
st_numscalar("r(sum)", colsum(x))
}
end

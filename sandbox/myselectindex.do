

mata: mata clear

mata

function selectindex2(bool_pattern){
	if(rows(bool_pattern)==1){
		n = cols(bool_pattern)
		A = (1..n)
		return(select(A,bool_pattern))
	} else {
		n = rows(bool_pattern)
		A = (1::n)
		return(select(A,bool_pattern))
	}
}

a = (1,2,0,8,0,11,2,7)
selectindex(a)
selectindex2(a)

selectindex(a')
selectindex2(a')

end

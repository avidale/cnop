version 12

mata

mata clear

class P{
	real scalar a
	real scalar b
}

class P scalar pmaker(su, di){
	class P scalar p
	p.a	= (su + di) / 2
	p.b	= (su - di) / 2
	return(p)
}



end

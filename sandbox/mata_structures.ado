version 12



mata: mata clear

mata:

struct abc{
	real scalar a,b,c
}

function main(){
	struct abc scalar x
	x = maker()
	x.a
	x.b
	x.c
}

struct abc scalar maker (){
	struct abc scalar x
	x.a=1
	x.b = 12
	return(x)
}

main()

end

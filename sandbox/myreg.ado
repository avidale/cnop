version 13

drop _all
capture program drop myreg 
capture program drop mypred

mata: mata clear

run C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\generate_CNOP.ado

mata:

n = 500; k = 5; ncat = 7; infcat = 4;
beta 	= (2 \ 3 \ 0 \ 0 \ 0.0);	 	a = (-0.4 \ 1.2);
gammap 	= (0 \ 1 \ 2 \ 0 \ 0.1);	mup = (-5 \ -2 \ 3 );
gamman 	= (0 \ 0 \ 3 \ 5 \ 0.2);	mun = (-4 \ -1 \ 3);
genCNOP(x=., y=., q=., n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun)
tmp = st_addvar("double", ("x":+strofreal(1..k), "y"))
st_addobs(n)
st_view(data = ., ., .)
data[,] = (x,y)

end

program mypred
	version 13
	syntax name [if] [in] 
	marksample touse
	local newVar = "`1'"
	mat b = e(b)
	local columnNames: colfullnames b
	tokenize `columnNames'
	gen `newVar' = b[1,1] + b[1,2] * `2'
end


program myreg, eclass
	version 13
	syntax varlist(min=2 max=2) [if] [in] 
	marksample touse
	
	matrix input b = (1.1, 2.3)
	matrix input V = (9, 1 \ 1, 4)
	
	matrix colnames b = _cons `2'
	matrix colnames V = _cons `2'
	matrix rownames V = _cons `2'

	ereturn post b V, esample(`touse')
	ereturn local predict "mypred"
	ereturn local cmd "myreg"
	ereturn display
end

myreg y x1

margins, dydx(x1)

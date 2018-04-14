drop _all
mata: mata clear

cd "C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox"


import delimited C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\test_cnop.csv, delimiter(";") 
replace y = y + 2
replace y2 = y2 + 3

run C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\stata_wrappers.ado

cnop y x1, infcat(2) zp(x2) zn(x2)

cnop y2 x1, infcat(3) zp(x2) zn(x2)

mata

st_view(y=., ., "y")
st_view(y2=., ., "y2")
st_view(x1=., ., "x1")
st_view(x2=., ., "x2")

st_view(lnl0=., ., "lnl0")


//lol = MLcnop(params, x, zp, zn, q, ncat, infcat) 

params = (0,0,1,0,0,0,0)
q = (y:==1), (y:==2), (y:==3)
x = x1
zp = x2
zn = x2
ncat = 3
infcat = 2
lol = MLcnop(params, x, zp, zn, q, ncat, infcat)
sum(lol)
lol-lnl0, y




n	= rows(x)
	kx	= cols(x)
	kzp	= cols(zp)
	kzn	= cols(zn)
ncatp = ncat - infcat
	ncatn = infcat - 1
	
	_cnop_params(params, kx, kzp, kzn, ncatp, ncatn, b = ., a = ., gp = ., mup = ., gn = ., mun = .)


xb = x*b
	zgp = zp*gp
	zgn = zn*gn
	
	p1n = normal(a[1]:-xb)
	p1p = 1 :- normal(a[2]:-xb)
	p10 = 1 :- p1p :- p1n


	p2p = normal(mup'[J(n,1,1),] :- zgp) , J(n,1,1)
	p2p[ ,2::(ncatp+1)] = p2p[,2::(ncatp+1)] - p2p[,1::ncatp]
	
	p2n = normal(mun'[J(n,1,1),] :- zgn) , J(n,1,1)
	p2n[,2::(ncatn+1)] = p2n[,2::(ncatn+1)] - p2n[,1::ncatn]


params = (0.4919,0.7492,2.5205,-0.8789,-0.1626	,2.6263,3.2775)


b^	0.5556			a^	1.8013	2.5011
g+^		-0.9698		mu+^	-0.2368	0.6061
g-^		0.7411		mu-^	-0.4636	0.6631


q2 = (y2:==1), (y2:==2), (y2:==3), (y2:==4), (y2:==5)
params = (0.5556,1.8013,2.5011,-0.9698,-0.2368, 0.6061, 0.7411, -0.4636, 0.6631)
x = x1
zp = x2
zn = x2
ncat = 5
infcat = 3
lol = MLcnop(params, x, zp, zn, q2, ncat, infcat)
sum(lol)

	
end



use "C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\testing\stat_data.dta", clear

cnop y x1 x3, infcat(0) zp(x2 x3) zn(x2 x3)

version 12

mata
/**/
function genZIOP(x, y, q, n, k, ncat, infcat, beta, a, gamma, mu) {
	x	= rnormal(n,k,0,1)
	q = J(n, ncat, 0)
	y	= J(n, 1, 1)
	nonz = (x * beta + rnormal (n,1,0,1) :> a)
	xg	= x * gamma + rnormal (n,1,0,1)
	for( i = 1; i <= ncat-1; i++){
		y = y + (xg :> mu[i])
	}
	y = y:*nonz+ infcat:*(1:-nonz);
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== i)
	}
}
//
function genMIOPR(x, y, q, n, k, ncat, infcat, beta, a, gamma, mu) {
// INCOMPLETE: include Z vars
	x	= rnormal(n,k,0,1)
	q = J(n, ncat, 0)
	y	= J(n, 1, 1)
	nonz = (x * beta + rnormal (n,1,0,1) :> a)
	xg	= x * gamma + rnormal (n,1,0,1)
	for( i = 1; i <= ncat-1; i++){
		y = y + (xg :> mu[i])
	}
	y = y:*nonz+ infcat:*(1:-nonz);	
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== i)
	}
}
//

function genCNOP(x, y, q, n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun){
	// here x = z always
	x	= rnormal(n,k,0,1)
	q = J(n, ncat, 0)
	y	= J(n, 1, 1)
	int0 = x * beta + rnormal (n,1,0,1)
	intp = x * gammap + rnormal (n,1,0,1)
	intn = x * gamman + rnormal (n,1,0,1)
	rp = (int0 :> a[2])
	rn = (int0 :< a[1])
	r0 = 1:-rp:-rn
	y = y + (rp+r0):*(infcat-1)
	
	for( i = 1; i <= infcat-1; i++){
		y = y + (intn :> mun[i]):*rn
	}
	for( i = 1; i <= ncat-infcat; i++){
		y = y + (intp :> mup[i]):*rp
	}
	
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== i)
	}
}

function genCNOPC(x, y, q, n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun, rhon, rhop){
	// here x = z always
	x	= rnormal(n,k,0,1)
	q = J(n, ncat, 0)
	y	= J(n, 1, 1)
	err0 = rnormal (n,1,0,1)
	errp = err0 * rhop + rnormal (n,1,0,1) * sqrt(1 - rhop^2)
	errn = err0 * rhon + rnormal (n,1,0,1) * sqrt(1 - rhon^2)
	int0 = x * beta + err0
	intp = x * gammap + errp
	intn = x * gamman + errn
	rp = (int0 :> a[2])
	rn = (int0 :< a[1])
	r0 = 1:-rp:-rn
	y = y + (rp+r0):*(infcat-1)
	
	for( i = 1; i <= infcat-1; i++){
		y = y + (intn :> mun[i]):*rn
	}
	for( i = 1; i <= ncat-infcat; i++){
		y = y + (intp :> mup[i]):*rp
	}
	
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== i)
	}
}


function genNOP(x, y, q, n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun){
	x	= rnormal(n,k,0,1)
	q = J(n, ncat, 0)
	y	= J(n, 1, 1)
	int0 = x * beta + rnormal (n,1,0,1)
	intp = x * gammap + rnormal (n,1,0,1)
	intn = x * gamman + rnormal (n,1,0,1)
	rp = (int0 :> a[2])
	rn = (int0 :< a[1])
	r0 = 1:-rp:-rn
	y = y + (rp+r0):*(infcat-1)
	y = y + rp:*1
	
	for(i = 1; i < infcat-1; i++){
		y = y + (intn :> mun[i]):*rn
	}
	for(i = 1; i < ncat-infcat; i++){
		y = y + (intp :> mup[i]):*rp
	}
	
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== i)
	}
}

function genNOPC(x, y, q, n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun, rhon, rhop){
	x	= rnormal(n,k,0,1)
	q = J(n, ncat, 0)
	y	= J(n, 1, 1)
	err0 = rnormal (n,1,0,1)
	errp = err0 * rhop + rnormal (n,1,0,1) * sqrt(1 - rhop^2)
	errn = err0 * rhon + rnormal (n,1,0,1) * sqrt(1 - rhon^2)
	int0 = x * beta + err0
	intp = x * gammap + errp
	intn = x * gamman + errn
	rp = (int0 :> a[2])
	rn = (int0 :< a[1])
	r0 = 1:-rp:-rn
	y = y + (rp+r0):*(infcat-1)
	y = y + rp:*1
	
	for( i = 1; i < infcat-1; i++){
		y = y + (intn :> mun[i]):*rn
	}
	for( i = 1; i < ncat-infcat; i++){
		y = y + (intp :> mup[i]):*rp
	}
	
	for(i=1;i<=ncat; i++){
		q[.,i]=(y :== i)
	}
}





function simulateCNOPC(x, n, k, ncat, infcat, beta, a, gammap, mup, gamman, mun, rhon, rhop) {
	y	= J(n, 1, 1)
	err0 = rnormal (n,1,0,1)
	errp = err0 * rhop + rnormal (n,1,0,1) * sqrt(1 - rhop^2)
	errn = err0 * rhon + rnormal (n,1,0,1) * sqrt(1 - rhon^2)
	int0 = x * beta + err0
	intp = x * gammap + errp
	intn = x * gamman + errn
	rp = (int0 :> a[2])
	rn = (int0 :< a[1])
	r0 = 1:-rp:-rn
	y = y + (rp+r0):*(infcat-1)
	
	for( i = 1; i <= infcat-1; i++){
		y = y + (intn :> mun[i]):*rn
	}
	for( i = 1; i <= ncat-infcat; i++){
		y = y + (intp :> mup[i]):*rp
	}
	
	return(y)
}


end

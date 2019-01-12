

/***
Before execution, we need to specify:
	DGP
	n
	start_iter
	sim_iter
	quiet

For example:
	DGP	= "NOPC"
	n	= 200
	start_iter	= 1
	sim_iter	= 10*1000
	quiet	= 1

Before execution, we need to call:
	CNOPishModel_definition.ado
	gradients.ado
	inflatedOP_estimation_routines.ado
	simulation_routines.do
***/


mata


generate_params(DGP, param_true=.,  beta=., a=., gammap=., gamman=., mup=., mun=., ron=., rop=., gamma=., mu=., ro=.)

ready = 0
con = 0

for(it = start_iter; it <= 50000; it++){
	it, ready, con
	
	/* Generate data */
	it2 = 0
	rseed(42+it)
	
	x1	= rnormal(n,1,0,1) :+ 2
	x2	= rnormal(n,1,0,1)
	x3	= runiform(n,1)
	x3	= (x3 :> 0.70) :- (x3 :< 0.30)

	e0 = rnormal(n,1,0,1)
	e1 = rnormal(n,1,0,1)
	e2 = rnormal(n,1,0,1)
	
	x	= x1, x2
	zp	= x2, x3
	zn	= x1, x3

	ncat = 5
	infcat = 3
	
	one_simulation(DGP, y=., n, ncat, infcat, x, zp, zn, e0, e1, e2, beta, a, gammap, gamman, mup, mun, ron, rop, gamma, mu, ro)
	
	y_hist = sum(y:==1), sum(y:==2), sum(y:==3), sum(y:==4), sum(y:==5)
	
	if (min(y_hist/rows(y))<0.06) {
		"bad data generated, continue another y"
		continue
	}
	
	if (min(y_hist) > 0) {
		estimate_and_get_params_v2(DGP, p=., s=., me=., mese = ., pr = ., prse = ., conv = ., etime = ., eiter = ., y=y, x=x, zp=zp, zn=zn, infcat=infcat, quiet=quiet)
	} else {
		"bad y hist"
		get_null_params_v2(DGP, p=., s=., me=., mese = ., pr = ., prse = ., conv = ., etime = ., eiter = ., y=y, x=x, zp=zp, zn=zn, infcat=infcat, quiet=quiet)
	}
	if (conv == 1) {
		"another successful row added!"
		if (ready == 0) {
			allpa = p
			allse = s
			allme = me
			allms = mese
			allpr = pr
			allps = prse
			allco = conv
			allet = etime
			allyh = y_hist
			allei = eiter
			allit = it
		} else {
			allpa = allpa \ p
			allse = allse \ s
			allme = allme \ me
			allms = allms \ mese
			allpr = allpr \ pr
			allps = allps \ prse
			allco = allco \ conv
			allet = allet \ etime
			allyh = allyh \ y_hist
			allei = allei \ eiter
			allit = allit \ it
		}
		ready = rows(allit)
		con = sum(allco)
		if(ready >= sim_iter){
			break
		}
	} else {
		"not converged, continue another y"
	}
}

it
rows(sim_iter)

param_true,  mean(allpa)'
mean(allco)

sim = allpa, allse, allme, allms, allpr, allps, allco, allet, allyh, allei, allit
/*
sim = allpa, allse, allme, allms, allpr, allps, allco, allet, allyh
*/

/* */
fname = "mc_" + strlower(DGP) + "_" + strofreal(n) + "_partial_"+strofreal(start_iter)+"to"+strofreal(sim_iter)+".matamatrix"
fh = fopen(fname, "w")
fputmatrix(fh, sim)
fclose(fh)


end



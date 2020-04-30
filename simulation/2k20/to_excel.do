/***
Before execution, we need to specify:
	dgp
	n
	start_iter
	sim_iter

For example:
DGP	= "NOPC"
n	= 200
start_iter = 1
sim_iter = 1000

Before execution, we need to call:
	CNOPishModel_definition.ado
	gradients.ado
	inflatedOP_estimation_routines.ado
	simulation_routines.do
***/

mata

dgp=DGP

signature = strlower(DGP) + "_" + strofreal(n) + "_partial_"+strofreal(start_iter)+"to"+strofreal(sim_iter)

/* READ AND PARSE FILE */
fname = "2k20/results/mc_" + signature +".matamatrix"
fh = fopen(fname, "r")
sim = fgetmatrix(fh)
fclose(fh)


fname = "2k20/results/mc_boot_" + signature +".matamatrix"
fh = fopen(fname, "r")
sim_boot = fgetmatrix(fh)
fclose(fh)


generate_params(dgp, param_true=.,  beta=., a=., gammap=., gamman=., mup=., mun=., ron=., rop=., gamma=., mu=., ro=.)

parlen = rows(param_true)



allpa = sim[,(1          )..(parlen     )]
allse = sim[,(parlen  + 1)..(parlen*2   )]
allme = sim[,(parlen*2+01)..(parlen*2+15)]
allms = sim[,(parlen*2+16)..(parlen*2+30)]
allpr = sim[,(parlen*2+31)..(parlen*2+35)]
allps = sim[,(parlen*2+36)..(parlen*2+40)]
allco = sim[,(parlen*2+41)]
allet = sim[,(parlen*2+42)]
allyh = sim[,(parlen*2+43)..(parlen*2+47)]



/* GENERATE TRUE P AND ME */

get_true_me_p(dgp, param_true, true_me = ., true_pr = .)

/* FILTER CORRECT ITERATIONS */

minshare  = rowmin(allyh) :/ rowsum(allyh)
filter = minshare:>=0.05 :& allco
okrows = selectindex(filter)
rows(okrows)



/* PROCESS BOOTSTRAP ESTIMATES */

nboot = rows(sim_boot) / rows(sim)

/* for each big iteration, we calculate variance of bootstrap estimates */
bootvars = J(rows(sim), parlen+5+15, 999)
bootvars = J(rows(sim), parlen+5+15, 999)
for(i=0; i< rows(sim); i++) {
	slice = sim_boot[1+i*nboot..(i+1)*nboot,4..cols(sim_boot)]
	slice_var = colsum((slice :- colsum(slice):/nboot):^2) :/ (nboot-1)
	bootvars[i+1,] = slice_var
}
bootstds = bootvars:^0.5

/* sim_boot = ready, it, booti, boot_p, boot_me, boot_pr */
bootspa = bootstds[,(1          )..(parlen     )]
bootsme = bootstds[,(parlen+1   )..(parlen+15  )]
bootspr = bootstds[,(parlen+16  )..(parlen+20  )]



/* MAKE ALL COMPARISONS */

cv = invnormal(.975)


meanparams = mean(allpa[okrows,])'
rmse = (mean((allpa[okrows,]:-param_true'):^2)'):^0.5

meanse = mean(allse[okrows,])'
medianse = colMedians(allse[okrows,])'
realse = (mean((allpa[okrows,]:-meanparams'):^2)'):^0.5
realse_boot = mean(bootspa[okrows,])'


meancoverage = calc_coverage(param_true, allpa[okrows,], allse[okrows,], cv)
boot_meancoverage = calc_coverage(param_true, allpa[okrows,], bootspa[okrows,], cv)



p_meanparams = mean(allpr[okrows,])'
p_rmse = (mean((allpr[okrows,]:-true_pr):^2)'):^0.5

p_meanse = mean(allps[okrows,])'
p_medianse = colMedians(allps[okrows,])'
p_realse = (mean((allpr[okrows,]:-p_meanparams'):^2)'):^0.5
p_realse_boot = mean(bootspr[okrows,])'

p_cil = allpr[okrows,] - cv:*allps[okrows,]
p_ciu = allpr[okrows,] + cv:*allps[okrows,]

p_coverage = (true_pr' :> p_cil') :* (true_pr' :< p_ciu')
p_meancoverage = mean(p_coverage')'
p_meancoverage = calc_coverage(true_pr', allpr[okrows,], allps[okrows,], cv)
p_boot_meancoverage = calc_coverage(true_pr', allpr[okrows,], bootspr[okrows,], cv)


m_meanparams = mean(allme[okrows,])'
m_rmse = (mean((allme[okrows,]:-true_me):^2)'):^0.5

m_meanse = mean(allms[okrows,])'
m_medianse = colMedians(allms[okrows,])'
m_realse = (mean((allme[okrows,]:-m_meanparams'):^2)'):^0.5
m_realse_boot = mean(bootsme[okrows,])'

m_meancoverage = calc_coverage(true_me', allme[okrows,], allms[okrows,], cv)
m_boot_meancoverage = calc_coverage(true_me', allme[okrows,], bootsme[okrows,], cv)




/*
todo: 
1) estimate standard error by 

- Coverage rate
- Bias of standard error estimates
*/




/* WRITE COMPARISON TO EXCEL */

excel = xl()

if (fileexists("2k20/results/sim_results.xlsx")) {
    excel.load_book("2k20/results/sim_results.xlsx")
} else {
	excel.create_book("2k20/results/sim_results", "test", "xlsx")
}

excel.set_mode("open")

sheetname = signature
if(max(excel.get_sheets():==sheetname)==0){
	excel.add_sheet(sheetname)
}
excel.set_sheet(sheetname)
/*
/* decided not to cleer, because additional formulas may be here */
excel.clear_sheet(sheetname)
*/

erow = 1

excel.put_string(erow, 1, "ESTIMATION RESULTS for " + sheetname)
erow = erow + 2


excel.put_string(erow, 1, "Converged: ")
excel.put_number(erow, 2, sum(allco))
erow=erow+1
excel.put_string(erow, 1, "Filtered:")
excel.put_number(erow, 2, sum(filter))
erow = erow + 2


excel.put_string(erow, 1, "PARAMETERS " + sheetname)
erow = erow + 1
excel.put_string(erow, 1, ("TRUE", "mean", "mean se", "real se", "real2mean se", "real2median se", "rmse", "coverage", "boot_cov", "boot_se", "real2boot_se"))
erow=erow+1
excel.put_number(erow, 1, (param_true,meanparams,meanse,realse,(realse :/ meanse),(realse :/ medianse), rmse,meancoverage, boot_meancoverage, realse_boot, (realse :/ realse_boot) ) )
erow = erow + rows(param_true) + 1

excel.put_string(erow, 1, "PROBABILITIES " + sheetname)
erow = erow + 1
excel.put_string(erow, 1, ("TRUE", "mean", "mean se", "real se", "real2mean se", "real2median se", "rmse", "coverage", "boot_cov", "boot_se", "real2boot_se"))
erow=erow+1
excel.put_number(erow, 1, (true_pr',p_meanparams,p_meanse,p_realse,(p_realse :/ p_meanse),(p_realse :/ p_medianse), p_rmse,p_meancoverage, p_boot_meancoverage, p_realse_boot, (p_realse :/ p_realse_boot) ) )
erow = erow + rows(true_pr') + 1

excel.put_string(erow, 1, "MARGINAL EFFECTS " + sheetname)
erow = erow + 1
excel.put_string(erow, 1, ("TRUE", "mean", "mean se", "real se", "real2mean se", "real2median se", "rmse", "coverage", "boot_cov", "boot_se", "real2boot_se"))
erow=erow+1
excel.put_number(erow, 1, (true_me',m_meanparams,m_meanse,m_realse,(m_realse :/ m_meanse),(m_realse :/ m_medianse), m_rmse,m_meancoverage, m_boot_meancoverage, m_realse_boot, (m_realse :/ m_realse_boot) ) )
erow = erow + rows(true_me') + 1



excel.close_book()



end



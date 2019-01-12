/***
Before execution, we need to specify:
	DGP
	n

For example:
DGP	= "NOPC"
n	= 200

Before execution, we need to call:
	CNOPishModel_definition.ado
	gradients.ado
	inflatedOP_estimation_routines.ado
	simulation_routines.do
***/




mata

/* READ AND PARSE FILE */

fname = "mc_" + strlower(dgp) + "_" + strofreal(n) + "_partial_1to100.matamatrix"

fh = fopen(fname, "r")
sim = fgetmatrix(fh)
fclose(fh)


/* dirty patches */ 

if(dgp == "NOP" & n == 200){
	sim = select(sim, sim[,7]:>-8000)
}
if(dgp == "NOPC" & n == 200){
	sim = select(sim, sim[,6]:>-1000)
}
if(dgp == "CNOPC" & n == 500){
	sim = select(sim, sim[,12]:<1e11)
}
if(dgp == "MIOPR" & n == 200){
	sim = select(sim, abs(sim[,1]):<200)
}
if(dgp == "CNOP" & n == 1000){
	sim = select(sim, sim[,12*2+16+2]:<9)
}
if(dgp == "CNOP" & n == 500){
	sim = select(sim, sim[,12*2+16+2]:<15)
}



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

/* MAKE ALL COMPARISONS */

cv = invnormal(.975)


meanparams = mean(allpa[okrows,])'
rmse = (mean((allpa[okrows,]:-param_true'):^2)'):^0.5

meanse = mean(allse[okrows,])'
medianse = colMedians(allse[okrows,])'
realse = (mean((allpa[okrows,]:-meanparams'):^2)'):^0.5

cil = allpa[okrows,] - cv:*allse[okrows,]
ciu = allpa[okrows,] + cv:*allse[okrows,]

coverage = (param_true :> cil') :* (param_true :< ciu')
meancoverage = mean(coverage')'


p_meanparams = mean(allpr[okrows,])'
p_rmse = (mean((allpr[okrows,]:-true_pr):^2)'):^0.5

p_meanse = mean(allps[okrows,])'
p_medianse = colMedians(allps[okrows,])'
p_realse = (mean((allpr[okrows,]:-p_meanparams'):^2)'):^0.5

p_cil = allpr[okrows,] - cv:*allps[okrows,]
p_ciu = allpr[okrows,] + cv:*allps[okrows,]

p_coverage = (true_pr' :> p_cil') :* (true_pr' :< p_ciu')
p_meancoverage = mean(p_coverage')'


m_meanparams = mean(allme[okrows,])'
m_rmse = (mean((allme[okrows,]:-true_me):^2)'):^0.5

m_meanse = mean(allms[okrows,])'
m_medianse = colMedians(allms[okrows,])'
m_realse = (mean((allme[okrows,]:-m_meanparams'):^2)'):^0.5

m_cil = allme[okrows,] - cv:*allms[okrows,]
m_ciu = allme[okrows,] + cv:*allms[okrows,]

m_coverage = (true_me' :> m_cil') :* (true_me' :< m_ciu')
m_meancoverage = mean(m_coverage')'




/* WRITE COMPARISON TO EXCEL */

excel = xl()
/*
excel.create_book("sim_results", "test") 
*/

excel.load_book("sim_results")
excel.set_mode("open")

sheetname = strupper(dgp) + "_" + strofreal(n)
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
excel.put_string(erow, 1, ("TRUE", "mean", "mean se", "real se", "real2mean se", "real2median se", "rmse", "coverage"))
erow=erow+1
excel.put_number(erow, 1, (param_true,meanparams,meanse,realse,(realse :/ meanse),(realse :/ medianse), rmse,meancoverage) )
erow = erow + rows(param_true) + 1

excel.put_string(erow, 1, "PROBABILITIES " + sheetname)
erow = erow + 1
excel.put_string(erow, 1, ("TRUE", "mean", "mean se", "real se", "real2mean se", "real2median se", "rmse", "coverage"))
erow=erow+1
excel.put_number(erow, 1, (true_pr',p_meanparams,p_meanse,p_realse,(p_realse :/ p_meanse),(p_realse :/ p_medianse), p_rmse,p_meancoverage) )
erow = erow + rows(true_pr') + 1

excel.put_string(erow, 1, "MARGINAL EFFECTS " + sheetname)
erow = erow + 1
excel.put_string(erow, 1, ("TRUE", "mean", "mean se", "real se", "real2mean se", "real2median se", "rmse", "coverage"))
erow=erow+1
excel.put_number(erow, 1, (true_me',m_meanparams,m_meanse,m_realse,(m_realse :/ m_meanse),(m_realse :/ m_medianse), m_rmse,m_meancoverage) )
erow = erow + rows(true_me') + 1



excel.close_book()



end



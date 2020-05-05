mata

/*
These parameters should be specified before running the script:
*/

DGP = "NOP"
n = 200

start_iter_1 = 1
start_iter_2 = 31
sim_iter_1 = 30
sim_iter_2 = 20

/* creating the filenames from the parameters */

name1 = strlower(DGP) + "_" + strofreal(n) + "_partial_"+strofreal(start_iter_1)+"to"+strofreal(sim_iter_1)
name2 = strlower(DGP) + "_" + strofreal(n) + "_partial_"+strofreal(start_iter_2)+"to"+strofreal(sim_iter_2)

name3 = strlower(DGP) + "_" + strofreal(n) + "_partial_"+strofreal(start_iter_1)+"to"+strofreal(sim_iter_1 + sim_iter_2)

/*
exporting variable names for to_excel.do
*/
dgp = DGP
start_iter = start_iter_1
sim_iter = sim_iter_1 + sim_iter_2


/* READ AND PARSE FILE 1 */
fh = fopen("2k20/results/mc_" + name1 +".matamatrix", "r")
sim1 = fgetmatrix(fh)
fclose(fh)

fh = fopen("2k20/results/mc_boot_" + name1 +".matamatrix", "r")
sim_boot1 = fgetmatrix(fh)
fclose(fh)

/* READ AND PARSE FILE 2 */
fh = fopen("2k20/results/mc_" + name2 +".matamatrix", "r")
sim2 = fgetmatrix(fh)
fclose(fh)

fh = fopen("2k20/results/mc_boot_" + name2 +".matamatrix", "r")
sim_boot2 = fgetmatrix(fh)
fclose(fh)


/* merge ans save the results */

fh = fopen("2k20/results/mc_" + name3 +".matamatrix", "w")
fputmatrix(fh, sim1 \ sim2)
fclose(fh)


fh = fopen("2k20/results/mc_boot_" + name3 +".matamatrix", "w")
fputmatrix(fh, sim_boot1 \ sim_boot2)
fclose(fh)

end

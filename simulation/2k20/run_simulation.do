mata: mata clear

cd "C:\Users\david\documents\cnop\sandbox"

set more off

run myselectindex.do
run CNOPishModel_definition.ado
run gradients.ado
run inflatedOP_estimation_routines.ado

cd "C:\Users\david\documents\cnop\simulation"
run 2k20/simulation_routines.do



mata

DGP	= "CNOPC"
n	= 200
n_boot = 7
start_iter	= 1
sim_iter	= 9
quiet	= 1

min_y_pct = 0.06
min_boot_y_pct = 0.06

end

do 2k20/simulation_script_v6.do


mata


end








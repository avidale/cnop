mata: mata clear

cd "C:\Users\david\yandexdisk\hsework\gauss-mata\cnop\sandbox"

set more off

run myselectindex.do
run CNOPishModel_definition.ado
run gradients.ado
run inflatedOP_estimation_routines.ado

cd "C:\Users\david\yandexdisk\hsework\gauss-mata\cnop\simulation"
run 2k20/simulation_routines.do



mata

DGP	= "NOP"
n	= 500
n_boot = 10
start_iter	= 1
sim_iter	= 15
quiet	= 1

min_y_pct = 0.06
min_boot_y_pct = 0.06

end

do 2k20/simulation_script_v6.do

do 2k20/to_excel.do

mata


end








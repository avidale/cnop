mata: mata clear

cd "C:\Users\ddale\YandexDisk\hsework\gauss-mata\cnop\sandbox"

set more off


run myselectindex.do

run CNOPishModel_definition.ado

run gradients.ado

run inflatedOP_estimation_routines.ado
/**/

cd "C:\Users\ddale\YandexDisk\hsework\gauss-mata\cnop\simulation"
run simulation_routines.do



mata

DGP	= "NOP"
n	= 200
start_iter	= 1
sim_iter	= 100
quiet	= 1

end

do simulation_script_v5.do


mata


end








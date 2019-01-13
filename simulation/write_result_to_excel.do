mata: mata clear

cd "C:\Users\ddale\YandexDisk\hsework\gauss-mata\cnop\sandbox"

run myselectindex.do
run CNOPishModel_definition.ado
run gradients.ado
run inflatedOP_estimation_routines.ado

cd "C:\Users\ddale\YandexDisk\hsework\gauss-mata\cnop\simulation"
run simulation_routines.do


mata

dgp			= "CNOP"
mdl			= "CNOP"
n			= 210
end

do simulation_comparison_script.do

mata: mata clear

cd "C:\Users\ddale\YandexDisk\hsework\gauss-mata\cnop\sandbox"

set more off

run myselectindex.do
run CNOPishModel_definition.ado	
run gradients.ado
run inflatedOP_estimation_routines.ado

cd "C:\Users\ddale\YandexDisk\hsework\gauss-mata\cnop\application"

use rate_change.dta, clear

cd "C:\Users\ddale\YandexDisk\hsework\gauss-mata\cnop\simulation"
run simulation_routines.do

mata

DGP	= "CNOP"
MDLS = "OP", "CNOP"
start_iter	= 1
sim_iter	= 3
times_per_iter = 3
bootstrap = 1
quiet	= 1
MIN_CLASS_PERCENTAGE = 0.02
MIN_CLASS_COUNT = 8
repeat_dataset = 1

fname = "sim_compare_results_2020.xlsx"
sheetname = "MC results (20200123)"

end

do simulation_script_different_dgp.do
do bisimulation_to_excel2.do

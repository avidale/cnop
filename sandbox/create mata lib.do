

cd "C:\Users\David\YandexDisk\hsework\gauss-mata\cnop\sandbox\"

// todo: make sure the library is not present in BASE

// ! del "C:\Program Files (x86)\Stata13\ado\base\l\lziop.mlib"

mata: mata clear

// define CNOPModel structure
run CNOPishModel_definition.ado
// define all auxiliary functions for optimization and ME
run gradients.ado
// define functions for model estimation
run inflatedOP_estimation_routines.ado


mata: mata mlib create lziop, replace
mata: mata mlib add lziop *()
mata: mata describe using lziop
mata: mata clear

// ! move C:\Users\David\YandexDisk\hsework\gauss-mata\cnop\package\lziop.mlib C:\Program Files (x86)\Stata13\ado\base\l\

// mata: mata mlib index

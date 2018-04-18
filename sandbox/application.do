
// clear workspace and load the needed functions
drop _all
mata: mata clear
// make the current directory working directory (change the address to that on your own computer!)
//cd C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox\
cd "C:\Users\user\Documents\Dale\Our paper on Github\cnop\sandbox"
// load the required procedures
run stata_wrappers.ado


import delimited Data_for_application.csv, clear 

cnop y5 pb_l pb_t spread houst gdp in 3/193, zn(spread gdp) zp(pb_t spread gdp) infcat(3)
cnop y3 pb spread houst gdp in 3/193, zn(spread gdp) zp(pb spread) infcat(2) correlated
nop y5 pb spread houst gdp in 3/193, zn(spread gdp) zp(pb spread) infcat(3)
nop y5 pb spread houst gdp in 3/193, zn(spread gdp) zp(pb spread) infcat(3) correlated

gen pb2 = pb
replace pb2 = 1 if pb != 0

miop y5 pb2 spread houst gdp in 3/193, z(pb_t pb_l spread houst gdp) infcat(3) correlated



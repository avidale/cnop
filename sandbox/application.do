
// clear workspace and load the needed functions
drop _all
mata: mata clear
// make the current directory working directory (change the address to that on your own computer!)
//cd C:\Users\David\YandexDisk\hsework\gauss-mata\cnop\sandbox\
cd "C:\Users\user\Documents\Dale\Our paper on Github\cnop\sandbox"
// load the required procedures
run stata_wrappers.ado


import delimited Data_for_application.csv, clear 

gen spreada = abs(spread)
set more off
nop y5 pb spread houst gdp in 5/214, xn(spread ) xp(pb spread) infcat(3)
set more off
nop y5 pb spread houst gdp in 5/214, xn(spread gdp) xp(pb spread) infcat(3) endoswitch
ziop2 y3 spreada  in 5/214, x(pb houst gdp ) infcat(2)
ziop2 y3 spreada  in 5/214, x(pb houst gdp ) infcat(2) endoswitch
ziop3 y5 pb spread houst gdp in 5/214, xn(spread gdp) xp(pb spread) infcat(3)
ziop3 y5 pb spread houst gdp in 5/214, xn(spread gdp) xp(pb spread) infcat(3) endoswitch

set more off
ziop2 y3 spreada in 5/214, x(houst pb gdp ) infcat(2) endoswitch



predict yfit
predict pr0, zeros
predict pr, regime
predict emode, output(mode)
predict emean, output(mean)
predict pcum, output(cum)

ziopmargins
ziopmargins, zeros
ziopmargins, regime

ziopprobabilities
ziopprobabilities, zeros
ziopprobabilities, regime

ziopcontrasts, at(pb=1) to(pb=0)
ziopcontrasts, at(pb=1) to(pb=0) zeros
ziopcontrasts, at(pb=1) to(pb=0) regime

ziopconfusion

//ziop3 y5 pb_l pb_t spread houst gdp in 5/214, xn(spread gdp) xp(pb_t spread gdp) infcat(3)
//ziop3 y5 pb_l pb_t spread houst gdp in 5/214, xn(spread gdp) xp(pb_t spread) infcat(3) endoswitch

gen yz = y5 -3
replace yz = 0 if y5 < 3
gen pb2 = pb
replace pb2 = 1 if pb != 0

gen spread_p = spread
replace spread_p = 0 if spread < 0
gen spread_n = spread
replace spread_n = 0 if spread > 0
gen y_L = y5_01
replace y_L = 1 if y5_01 != 3
replace y_L = 0 if y5_01 == 3
gen y_Lp = 1
gen y_Ln = 1
replace y_Lp = 0 if y5_01 < 4
replace y_Ln = 0 if y5_01 > 2
gen houst_g = houst - 1.4664
gen houst_p = houst_g
gen houst_n = houst_g
replace houst_p = 0 if houst_g < 0
replace houst_n = 0 if houst_g > 0


//ziop2 y5 pb_l pb_t in 5/214, x(pb spread houst gdp) infcat(3) endoswitch


// vuong example
ziop3 y5 pb spread houst gdp in 5/214, xn(spread gdp) xp(pb spread) infcat(3)
set more off
est store firstmodel
nop y5 pb spread houst gdp in 5/214, xn(spread gdp) xp(pb spread) infcat(3)
set more off
est store secondmodel
ziopvuong firstmodel secondmodel

//confusion example
ziop3 y5 pb spread houst gdp in 5/214, xn(spread gdp) xp(pb spread) infcat(3)
set more off
ziopclassification

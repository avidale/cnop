
// clear workspace and load the needed functions
drop _all
mata: mata clear
// make the current directory working directory (change the address to that on your own computer!)
//cd C:\Users\David\YandexDisk\hsework\gauss-mata\cnop\sandbox\
cd "C:\Users\user\Documents\Dale\Our paper on Github\cnop\sandbox"
// load the required procedures
run stata_wrappers.ado

import delimited Data_for_application.csv, clear 

set more off
nop rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0)
set more off
nop rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) endoswitch
set more off
ziop2 rate_change spread_u spread_d houst_u houst_d, x(spread pb houst gdp ) infcat(0)
set more off
ziop2 rate_change spread_u spread_d  gdp_u gdp_d houst_u houst_d, x(spread pb houst gdp ) infcat(0) endoswitch
set more off
ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0)
set more off
ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) endoswitch

set more off
ziop3 or_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) endoswitch

set more off
quietly ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) endoswitch
predict proby

predict przeros, zeros
predict prregim, regime
predict emode, output(mode)
predict emean, output(mean)
predict pcum, output(cum)

quietly ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) endoswitch
ziopmargins
ziopmargins, zeros
ziopmargins, regime

ziopprobabilities, at (pb=1, spread=0.426, houst=1.6, gdp=6.8)



ziopprobabilities, zeros
ziopprobabilities, regime

ziopcontrasts, at(pb=1) to(pb=0)
ziopcontrasts, at(pb=1) to(pb=0) zeros
ziopcontrasts, at(pb=1) to(pb=0) regime

// vuong example
quietly ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(0) endoswitch
est store firstmodel
quietly ziop2 rate_change houst_u houst_d spread_u spread_d, x( houst spread pb  gdp ) infcat(0) endoswitch
est store secondmodel
ziopvuong firstmodel secondmodel

//classification example
set more off
ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(0) endoswitch
ziopclassification


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
nop rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(3)
set more off
nop rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(3) endoswitch
set more off
ziop2 rate_change spread_u spread_d houst_u houst_d, x(spread pb houst gdp ) infcat(3)
set more off
ziop2 rate_change spread_u spread_d  gdp_u gdp_d houst_u houst_d, x(spread pb houst gdp ) infcat(3) endoswitch
set more off
ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(3)
set more off
ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(3) endoswitch

set more off
ziop2 rate_change spread_u spread_d pb_u pb_d houst_u houst_d gdp_u gdp_d in 1/241, x(spread pb houst gdp ) infcat(3) endoswitch


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

// vuong example
set more off
ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(3) endoswitch
est store firstmodel
set more off
ziop2 rate_change houst_u houst_d spread_u spread_d, x( houst spread pb  gdp ) infcat(3) endoswitch
est store secondmodel
ziopvuong firstmodel secondmodel

//classification example
set more off
ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(3) endoswitch
ziopclassification

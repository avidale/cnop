
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
oprobit rate_change spread pb houst gdp, nolog
estat ic

set more off
nop rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) nolog

set more off
nop rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) endoswitch nolog

set more off
ziop2 rate_change spread_u spread_d houst_u houst_d, x(spread pb houst gdp ) infcat(0)

set more off
ziop2 rate_change spread_u spread_d houst_u houst_d, x(spread pb houst gdp ) infcat(0) endoswitch

set more off
ziop2 rate_change spread pb houst gdp, x(spread pb houst gdp ) infcat(0) nolog

set more off
ziop2 rate_change spread_u spread_d  gdp_u gdp_d houst_u houst_d, x(spread pb houst gdp ) infcat(0) endoswitch nolog

set more off
ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) nolog

set more off
ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) endoswitch nolog


set more off
quietly ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) nolog
predict proby
predict proby, at(pb=1, spread=0.426, houst=1.6, gdp=6.8)

predict przeros, zeros
predict prregim, regime
predict emode, output(mode)
predict emean, output(mean)
predict pcum, output(cum)

predict pcum, output(cum) at(pb=1, spread=0.426, houst=1.6, gdp=6.8)

quietly ziop3 rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0) nolog
ziopmargins, at (pb=1, spread=0.426, houst=1.6, gdp=6.8)
ziopmargins
ziopmargins, zeros
ziopmargins, regime

ziopprobabilities, at (pb=1, spread=0.426, houst=1.6, gdp=6.8)
ziopprobabilities
ziopprobabilities, zeros at (pb=1, spread=0.426, houst=1.6, gdp=6.8)
ziopprobabilities, regime at (pb=1, spread=0.426, houst=1.6, gdp=6.8)

ziopcontrasts, at(pb=1, spread=0.426, houst=1.6, gdp=6.8) to(pb=0, spread=-1.394, houst=1.2, gdp=1.9)
ziopcontrasts, at(pb=1, spread=0.426, houst=1.6, gdp=6.8)
ziopcontrasts, at(pb=1) to(pb=0) zeros
ziopcontrasts, at(pb=1) to(pb=0) regime

// vuong example ZIOP-3 vs ZIOP-2
set more off
quietly ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(0)
est store ziop3model
set more off
ziop2 rate_change spread_u spread_d houst_u houst_d, x(spread pb houst gdp ) infcat(0)
est store ziop2model
ziopvuong ziop3model ziop2model


quietly ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(0) 
est store ziop3model
quietly ziop2 rate_change spread_u spread_d houst_u houst_d, x(spread pb houst gdp ) infcat(0) 
est store ziop2model
ziopvuong ziop3model ziop2model

quietly ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(0)
est store ziop3model
quietly oprobit rate_change spread pb houst gdp, nolog
est store opmodel
ziopvuong ziop3model opmodel

//classification example
set more off
ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(0)
ziopclassification

set more off
quietly nop rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0)
set more off
ziopclassification

// view help
view ../package/ziop.sthlp
view ../package/ziop_postestimation.sthlp


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
nop rate_change spread pb houst gdp, neg_indepvars(spread gdp) pos_indepvars(spread pb) infcat(0) nolog

set more off
nop rate_change spread pb houst gdp, neg_indepvars(spread gdp) pos_indepvars(spread pb) infcat(0) endo nolog


set more off
ziop2 rate_change spread pb houst gdp, indepvars(spread pb houst gdp ) infcat(0) nolog


set more off
ziop3 rate_change spread pb houst gdp, neg_indepvars(spread gdp) pos_indepvars(spread pb) infcat(0) nolog

set more off
ziop3 rate_change spread pb houst gdp, neg_indepvars(spread gdp) pos_indepvars(spread pb) infcat(0) endos nolog


set more off
quietly ziop3 rate_change spread pb houst gdp, neg_indepvars(spread gdp) pos_indepvars(spread pb) infcat(0) nolog

predict proby

predict p_zero, zeros
predict p_reg, regimes
tabstat p_zero* p_reg*, stat(mean)


predict choice, output(choice)
predict emean, output(mean)
predict pcum, output(cum)

predict pcum, output(cum) at(pb=1, spread=0.426, houst=1.6, gdp=6.8)

quietly ziop3 rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) infcat(0) nolog
ziopmargins, at (pb=1, spread=0.426, houst=1.6, gdp=6.8)
ziopmargins
ziopmargins, zeros
ziopmargins, regimes

ziopprobabilities, at (pb=1, spread=0.426, houst=1.6, gdp=6.8)
ziopprobabilities
ziopprobabilities, zeros at (pb=1, spread=0.426, houst=1.6, gdp=6.8)
ziopprobabilities, regimes at (pb=1, spread=0.426, houst=1.6, gdp=6.8)

ziopcontrasts, at(pb=1, spread=0.426, houst=1.6, gdp=6.8) to(pb=0, spread=-1.394, houst=1.2, gdp=1.9)
ziopcontrasts, at(pb=1, spread=0.426, houst=1.6, gdp=6.8)
ziopcontrasts, at(pb=1) to(pb=0) zeros
ziopcontrasts, at(pb=1) to(pb=0) regime

// vuong example ZIOP-3 vs ZIOP-2
quietly ziop3 rate_change pb spread houst gdp, neg(spread gdp )pos(pb spread) infcat(0)
est store ziop3_model
quietly ziop2 rate_change spread pb houst gdp, ind(spread pb houst gdp ) infcat(0)
est store ziop2_model
ziopvuong ziop3_model ziop2_model


quietly ziop3 rate_change pb spread houst gdp, neg(spread gdp )pos(pb spread) infcat(0)
est store ziop3model
quietly oprobit rate_change spread pb houst gdp, nolog
est store opmodel
ziopvuong ziop3model opmodel

//classification example
<<<<<<< HEAD
quietly ziop3 rate_change pb spread houst gdp, xn(spread gdp )xp(pb spread) infcat(0)
ziopclassification

quietly ziop2 rate_change spread pb houst gdp, x(spread pb houst gdp ) infcat(0) nolog
ziopclassification

quietly nop rate_change spread pb houst gdp, xn(spread gdp) xp(spread pb) infcat(0)
=======
set more off
ziop3 rate_change pb spread houst gdp, neg(spread gdp )pos(pb spread) infcat(0)
ziopclassification

set more off
quietly nop rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) infcat(0)
set more off
>>>>>>> 3e79c29e7d37355a8c05efa37a46216276e3713c
ziopclassification

// view help
view ../package/ziop.sthlp
view ../package/ziop_postestimation.sthlp


set more off
ziop2 rate_change spread_u spread_d houst_u houst_d, ind(spread pb houst gdp ) infcat(0)

set more off
ziop2 rate_change spread_u spread_d houst_u houst_d, ind(spread pb houst gdp ) infcat(0) endoswitch

set more off
ziop2 rate_change spread_u spread_d  gdp_u gdp_d houst_u houst_d, ind(spread pb houst gdp ) infcat(0) endoswitch nolog

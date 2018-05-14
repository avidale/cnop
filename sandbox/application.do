
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
nop rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) inf(0) endo nolog vuong

set more off
ziop2 rate_change spread pb houst gdp, indepvars(spread pb houst gdp ) infcat(0) nolog

set more off
ziop3 rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) inf(0) nolog vuong

predict p_zero, zeros
predict p_reg, regimes
tabstat p_zero* p_reg*, stat(mean)

ziopmargins, at (pb=1, spread=0.426, houst=1.6, gdp=6.8)

ziopprobabilities, at (pb=1, spread=0.426, houst=1.6, gdp=6.8)

ziopcontrasts, at(pb=1, spread=0.426, houst=1.6, gdp=6.8) ///
               to(pb=0, spread=0.426, houst=1.6, gdp=6.8)

// vuong example ZIOP-3 vs ZIOP-2
quietly ziop3 rate_change pb spread houst gdp, neg(spread gdp ) pos(pb spread) inf(0)
est store ziop3_model
quietly ziop2 rate_change spread pb houst gdp, indepvars(spread pb houst gdp) inf(0)
est store ziop2_model
ziopvuong ziop3_model ziop2_model


//classification example

set more off
quietly ziop3 rate_change pb spread houst gdp, neg(spread gdp ) pos(pb spread) inf(0)
ziopclassification

set more off
quietly ziop2 rate_change spread pb houst gdp, indepvars(spread pb houst gdp ) infcat(0) nolog
ziopclassification

set more off
quietly nop rate_change spread pb houst gdp, neg_indepvars(spread gdp) pos_indepvars(spread pb) infcat(0)
ziopclassification


// view help
view ../package/ziop.sthlp
view ../package/ziop_postestimation.sthlp


{smcl}
{* *! version 0.0.1  08may2018}{...}
{title:Title}

{cmd:nop}   -- Nested ordered probit regression
{cmd:ziop2} -- Two-part zero-inflated ordered probit regression
{cmd:ziop3} -- Three-part zero-inflated ordered probit regression


{title:Syntax}

{cmd:ziop3} {depvar} {indepvars_reg} {ifin} {bind:[{cmd:,} {it:options}]}

{cmd:ziop2} {depvar} {indepvars} {ifin} {bind:[{cmd:,} {it:options}]}

{cmd:nop}   {depvar} {indepvars} {ifin} {bind:[{cmd:,} {it:options}]}


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}

{syntab :Second-stage covariates}
{synopt :{opth xp(varlist)}} list of covariates for positive response in {cmd:nop} and {cmd:ziop3} models; 
by default, it equals indepvars, the list of covariates for initial stage{p_end}
{synopt :{opth xn(varlist)}} list of covariates for negative response in {cmd:nop} and {cmd:ziop3} models; 
by default, it equals indepvars, the list of covariates for initial stage{p_end}
{synopt :{opth x(varlist)}} list of covariates for non-zero response in {cmd:ziop2} models; 
by default, it equals indepvars, the list of covariates for initial stage{p_end}

{syntab :Other model specification}
{synopt :{opt infcat(integer)}} value of the response variable that should be modeled as infated;
by default, it equals 0{p_end}
{synopt :{opt endoswitch}}  flag that errors in the first and second stages may be correlated,
forcing estimation of endogenous switching models{p_end}
{synopt :{opt robust}}  flag that variance-covariance estimator must be robust (based on
"sandwich") estimate{p_end}
{synopt :{opth cluster(varname)}}  clustering variable for robust variance estimator{p_end}
 
{syntab :Control of optimization}
{synopt :{opt initial(string)}}  whitespace-delimited list of initial parameter values for estimation,
in the following order: beta, alpha, gamma+, mu+, gamma-, mu-, rho+, rho-
{p_end}
{synopt :{opt nolog}} flag that intermediate results of optimization should not be displayed
{p_end}
{synoptline}

See {help ziop_postestimation:ziop postestimation} for features available after estimation.

{title:Description}

{cmd:nop} estimates a three-part nested ordered probit (NOP) regression of {depvar} on three possibly different sets of covariates: {indepvars_reg} in the regime equation, {pos_indepvars()} in the outcome equation conditional on the regime s=1, and {neg_indepvars()} in the outcome equation conditional on the regime s=-1.

An ordinal dependent variable depvar is assumed to take on at least five discrete ordinal values in the NOP model, at least two --- in the ZIOP-2 model, and at least three --- in the ZIOP-3 model. A list of the covariates in the regime equation indepvarsreg may be different from the lists of the covariates in the outcome equations.

{cmd:ziop3} estimates by ML the three-part cross-nested zero-inflated OP model with
possibly different sets of covariates in the regime and outcome equations and possibly 
endogenous switching among three latent regimes.

{cmd:ziop2} command estimates by ML the two-part cross-nested zero-inflated OP model with
possibly different sets of covariates in the regime and outcome equations and possibly 
endogenous switching among two latent regimes.



{title:Stored results}

{pstd}
{cmd:nop}, {cmd:ziop2} and {cmd:ziop3} store the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k_cat)}}number of categories{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(r2_p)}}pseudo-R-squared{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log likelihood, constant-only model{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:nop}, {cmd:ziop2} and {cmd:ziop3}, respectively{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(ll_obs)}}vector of observation-wise log-likelihood{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

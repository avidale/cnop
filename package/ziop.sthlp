{smcl}
{* *! version 0.0.1  08may2018}{...}
{title:Title}

{cmd:nop}   -- Nested ordered probit regression
{cmd:ziop2} -- Two-part zero-inflated ordered probit regression
{cmd:ziop3} -- Three-part zero-inflated ordered probit regression


{title:Syntax}

{cmd:ziop3} {depvar} {it:indepvars_reg} {ifin} {bind:[{cmd:,} {it:options}]}

{cmd:ziop2} {depvar} {it:indepvars_reg} {ifin} {bind:[{cmd:,} {it:options}]}

{cmd:nop}   {depvar} {it:indepvars_reg} {ifin} {bind:[{cmd:,} {it:options}]}


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}

{syntab :Model}
{synopt :{opth pos:_indepvars(varlist)}} independent variables in the outcome equation conditional on the regime s=1 for nonnegative outcomes in the {cmd:nop} and {cmd:ziop3} regressions; by default, it is identical to {it:indepvars_reg}, a list of the independent variables in the regime equation.{p_end}
{synopt :{opth neg:_indepvars(varlist)}} independent variables in the outcome equation conditional on the regime s=-1 for nonpositive outcomes in the {cmd:nop} and {cmd:ziop3} regressions; by default, it is identical to {it:indepvars_reg}, a list of the independent variables in the regime equation.{p_end}
{synopt :{opth ind:epvars(varlist)}} independent variables in the outcome equation of the {cmd:ziop2} regression; by default, it is identical to {it:indepvars_reg}, a list of the independent variables in the regime equation.{p_end}
{synopt :{opt inf:cat(choice)}} inflated choice -- value of the dependent variable in the regime s=0; by default, it equals 0.{p_end}
{synopt :{opt endo:switch}} use endogenous regime switching instead of default exogenous switching (regime switching is endogenous if the errors in the regime equation are correlated with the errors in the outcome equations, and exogenous otherwise).{p_end}

{syntab :SE/Robust}
{synopt :{opt robust}} use robust sandwich estimator of variance; the default estimator is based on the observed information matrix.{p_end}
{synopt :{opth cluster(varname)}} clustering variable for the clustered robust sandwich estimator of variance{p_end}

{syntab :Reporting}
{synopt :{opt vuong}} perform Vuong test against the conventional ordered probit model (not available for {cmd:ziop2}).{p_end}
 
{syntab :Maximization}
{synopt :{opt initial(string)}} whitespace-delimited list of the starting values of the parameters in the following order: gamma, mu, beta+, alpha+, beta-, alpha-, rho- and rho+ for the {cmd:nop} and {cmd:ziop3} regressions, and gamma, mu, beta, alpha and pho for the {cmd:ziop2} regression.{p_end}
{synopt :{opt nolog}} suppress the iteration log and intermediate results.{p_end}
{synoptline}

See {help ziop_postestimation:ziop postestimation} for features available after estimation.

{title:Description}

{cmd:nop} estimates a three-part nested ordered probit (NOP) model of ordinal variable {depvar}, which takes on at least five values, on three sets of independent variables: {it:indepvars_reg} in the regime equation, {cmd:pos_indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=1, and {cmd:neg_indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=-1.

{cmd:ziop2} estimates a two-part zero-inflated ordered probit (ZIOP-2) model of ordinal variable {depvar} on two sets of independent variables: {it:indepvars_reg} in the regime equation and {cmd:indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=1.

{cmd:ziop3} estimates a three-part zero-inflated ordered probit (ZIOP-3) model of ordinal variable {depvar}, which takes on at least three values, on three sets of independent variables: {it:indepvars_reg} in the regime equation, {cmd:pos_indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=1, and {cmd:neg_indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=-1.

The actual values taken on by the dependent variable are irrelevant, except that larger values are assumed to correspond to "higher" outcomes.

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

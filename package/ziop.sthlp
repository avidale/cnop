{smcl}
{* *! version 0.0.1  08may2018}{...}
{title:Title}

{pstd}{helpb nop}   -- Three-part nested ordered probit regression{p_end}
{pstd}{helpb ziop2} -- Two-part zero-inflated ordered probit regression{p_end}
{pstd}{helpb ziop3} -- Three-part zero-inflated ordered probit regression{p_end}


{title:Syntax}

{pstd}{cmd:ziop3} {depvar} {it:indepvars_reg} {ifin} {bind:[{cmd:,} {it:options}]}{p_end}

{pstd}{cmd:ziop2} {depvar} {it:indepvars_reg} {ifin} {bind:[{cmd:,} {it:options}]}{p_end}

{pstd}{cmd:nop}   {depvar} {it:indepvars_reg} {ifin} {bind:[{cmd:,} {it:options}]}{p_end}


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}

{syntab :Model}
{synopt :{opth pos:_indepvars(varlist)}} independent variables in the outcome equation conditional on the regime s=1 for nonnegative outcomes in the {cmd:nop} and {cmd:ziop3} regressions; by default, it is identical to {it:indepvars_reg}, a list of the independent variables in the regime equation.{p_end}
{synopt :{opth neg:_indepvars(varlist)}} independent variables in the outcome equation conditional on the regime s=-1 for nonpositive outcomes in the {cmd:nop} and {cmd:ziop3} regressions; by default, it is identical to {it:indepvars_reg}, a list of the independent variables in the regime equation.{p_end}
{synopt :{opth ind:epvars(varlist)}} independent variables in the outcome equation of the {cmd:ziop2} regression; by default, it is identical to {it:indepvars_reg}, a list of the independent variables in the regime equation.{p_end}
{synopt :{opt inf:cat(choice)}} value of the dependent variable in the regime s=0 (an inflated choice in {cmd:ziop2} and {cmd:ziop3} models; a neutral choice in {cmd:nop} model); by default, it equals 0.{p_end}
{synopt :{opt endo:switch}} use endogenous regime switching instead of default exogenous switching (regime switching is endogenous if the errors in the regime equation are correlated with the errors in the outcome equations, and exogenous otherwise).{p_end}

{syntab :SE/Robust}
{synopt :{opt robust}} use robust sandwich estimator of variance; the default estimator is based on the observed information matrix.{p_end}
{synopt :{opth cluster(varname)}} clustering variable for the clustered robust sandwich estimator of variance{p_end}

{syntab :Reporting}
{synopt :{opt vuong}} perform the Vuong test (Vuong 1989) against the conventional ordered probit (OP) model (not available for {cmd:ziop2}).{p_end}
 
{syntab :Maximization}
{synopt :{opt initial(string)}} whitespace-delimited list of the starting values of the parameters in the following order: gamma, mu, beta+, alpha+, beta-, alpha-, rho- and rho+ for the {cmd:nop} and {cmd:ziop3} regressions, and gamma, mu, beta, alpha and pho for the {cmd:ziop2} regression.{p_end}
{synopt :{opt nolog}} suppress the iteration log and intermediate results.{p_end}
{synoptline}

{pstd}See {help ziop_postestimation:ziop postestimation} for features available after estimation.{p_end}

{title:Description}

{pstd}{cmd:nop} estimates a three-part nested ordered probit (NOP) model of an ordinal variable {depvar}, which takes on at least five values, on three sets of independent variables: {it:indepvars_reg} in the regime equation, {cmd:pos_indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=1, and {cmd:neg_indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=-1 (Sirchenko 2013).{p_end}

{pstd}{cmd:ziop2} estimates a two-part zero-inflated ordered probit (ZIOP-2) model of an ordinal variable {depvar} on two sets of independent variables: {it:indepvars_reg} in the regime equation and {cmd:indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=1 (Harris and Zhao 2007; Brooks, Harris and Spencer 2012; Bagozzi and Mukherjee 2012).{p_end}

{pstd}{cmd:ziop3} estimates a three-part zero-inflated ordered probit (ZIOP-3) model of an ordinal variable {depvar}, which takes on at least three values, on three sets of independent variables: {it:indepvars_reg} in the regime equation, {cmd:pos_indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=1, and {cmd:neg_indepvars}{it:(varlist)} in the outcome equation conditional on the regime s=-1 (Sirchenko 2013).{p_end}

{pstd}The actual values taken on by the dependent variable are irrelevant, except that larger values are assumed to correspond to "higher" outcomes.{p_end}

{title:Examples}

{pstd}Setup
        . webuse rate_change

{pstd}Fit three-part nested ordered probit model with exogenous switching{p_end}
        . nop rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) inf(0)

{pstd}Fit three-part nested ordered probit model with endogenous switching and report Vuong test of NOP versus OP{p_end}
        . nop rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) inf(0) endo vuong

{pstd}Fit two-part zero-inflated ordered probit model with exogenous switching{p_end}
        . ziop2 rate_change spread pb houst gdp, ind(spread pb houst gdp) inf(0)

{pstd}Fit three-part zero-inflated ordered probit model with exogenous switching{p_end}
        . ziop3 rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) inf(0)

{pstd}Fit three-part zero-inflated ordered probit model with endogenous switching and report Vuong test of ZIOP-3 versus OP{p_end}
        . ziop3 rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) inf(0) endo vuong

{title:Stored results}

{pstd}
{cmd:nop}, {cmd:ziop2} and {cmd:ziop3} store the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
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

{title:References}

{pstd}Bagozzi, B. E., and B. Mukherjee. 2012. A mixture model for middle category inflation in ordered survey responses. {it:Political Analysis} 20: 369--386.{p_end}

{pstd}Brooks, R., M. N. Harris, and C. Spencer. 2012. Inflated ordered outcomes. {it:Economics Letters} 117: 683--686.{p_end}

{pstd}Harris, M. N., and X. Zhao. 2007. A zero-inflated ordered probit model, with an application to modelling tobacco consumption. {it:Journal of Econometrics} 141: 1073--1099.{p_end}

{pstd}Sirchenko, A. 2013. A model for ordinal responses with an application to policy interest rate. National Bank of Poland Working Paper No. 148.{p_end}

{pstd}Vuong, Q. H. 1989. Likelihood ratio tests for model selection and non-nested hypotheses. {it:Econometrica} 57: 307-333.{p_end}

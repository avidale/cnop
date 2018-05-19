{smcl}
{* *! version 0.0.1  08may2018}{...}
{title:Title}

{pstd}
{helpb ziop postestimation} -- Postestimation tools fo {cmd:nop}, {cmd:ziop2} and {cmd:ziop3}

{title:Postestimation commands}

{pstd}
The following postestimation commands are available after {cmd:nop}, {cmd:ziop2} and {cmd:ziop3}: 

{synoptset 20 notes}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
{synopt :{helpb ziop postestimation##predict:predict}}predicted probabilities and other predictions for all values of independent variables{p_end}
{synopt :{helpb ziop postestimation##ziopprobabilities:ziopprobabilities}}predicted probabilities for specified values of independent variables{p_end}
{synopt :{helpb ziop postestimation##ziopcontrasts:ziopcontrasts}}differences in predicted probabilities for specified values of independent variables{p_end}
{synopt :{helpb ziop postestimation##ziopmargins:ziopmargins}}marginal effects on probabilities for specified values of independent variables{p_end}
{synopt :{helpb ziop postestimation##ziopclassification:ziopclassification}}classification table and other goodness-of-fit measures{p_end}
{synopt :{helpb ziop postestimation##ziopvuong:ziopvuong}}Vuong test for non-nested hypotheses{p_end}
{synoptline}
{p2colreset}{...}


{marker predict}{...}
{title:Syntax for predict}

{pstd}
{cmd:predict} {varname} {ifin} [, {opt zeros} {opt regimes} {opt output(string)} ]

{synoptset 20 notes}{...}
{p2coldent :Option}Description{p_end}
{synoptline}
{synopt :{cmd:zeros}}indicates that the probabilities of different types of zeros (the outcomes in the inflated category specified in {opt infcat(choice)}), conditional on different regimes, must be predicted instead of the choice probabilities{p_end}
{synopt :{cmd:regimes}}indicates that the probabilities of the regimes must be predicted instead of the choice probabilities; this option is ignored if the option {cmd:zeros}} is used.{p_end}
{synopt :{opt output(string)}}specifies the different types of predictions. The possible options for {it:string} are: {it:choice} for reporting the predicted outcome (the choice with the largest predicted probability); {it:mean} for reporting the expected value of the dependent variable computed as a summation of i*Pr(y=i) across all choices i; and 
{it:cum} for predicting the cumulative choice probabilities such as Pr(y<=0), Pr(y<=1), ... . If {it:string} is not specified, the usual choice probabilities such as Pr(y=0), Pr(y=1), ... are predicted and saved into new variables with the {varname} prefix.{p_end}
{synoptline}
{p2colreset}{...}

{title:Description for predict}

{pstd}
    {cmd:predict} creates new variables containing predictions such as the predicted probabilities of the discrete choices, the regimes, the types of zeros conditional on the regime, the expected values of the dependent variable, the predicted choice (one with the largest predicted probability) for all sample values of the independent variables.

{title:Examples}

{pstd}Setup{p_end}
        . webuse rate_change
        . ziop3 rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) inf(0) endo

{pstd}Predicted probabilities of discrete choices{p_end}
        . predict pr_choice

{pstd}Predicted discrete choice (one with the largest probability){p_end}
        . predict pr_choice output(choice)

{pstd}Expected value of dependent variable{p_end}
        . predict pr_choice output(mean)

{pstd}Predicted cumulative probabilities of discrete choices{p_end}
        . predict pr_choice output(cum)

{pstd}Predicted probabilities of three types of zeros conditional on the regime{p_end}
        . predict pr_zero, zeros

{pstd}Predicted probabilities of three regimes{p_end}
        . predict pr_regime, regimes



{marker ziopprobabilities}{...}
{title:Syntax for ziopprobabilities}

{pstd}
{cmd:ziopprobabilities} [, {opt at(string)}  {opt zeros} {opt regimes} ]

{synoptset 20 notes}{...}
{p2coldent :Option}Description{p_end}
{synoptline}
{synopt :{cmd:zeros}}indicates that the probabilities of different types of zeros (the outcomes in the inflated category specified in {opt infcat(choice)}), conditional on different regimes, must be predicted instead of the choice probabilities.{p_end}
{synopt :{cmd:regimes}}indicates that the probabilities of the regimes must be predicted instead of the choice probabilities; this option is ignored if the option {cmd:zeros} is used.{p_end}
{synopt :{opt at(string)}}specifies for which values of the independent variables to estimate the predictions. If at() is used ({it:string} is a list of varname=value expressions, separated by commas), the predictions are estimated for these values and displayed without saving to the dataset. If some independent variable names are not specified, their median values are taken instead. If at() is not used, by default the predictions are estimated for the median values of the independent variables.{p_end}
{p2colreset}{...}

{title:Description for ziopprobabilities}

{pstd}
    {cmd:ziopprobabilities} shows the predicted probabilities estimated for the specified values of independent variables along with the standard errors.{p_end}

{title:Examples}

{pstd}Setup{p_end}
        . webuse rate_change
        . ziop3 rate_change spread pb houst gdp, neg(spread gdp) pos(spread pb) inf(0) endo

{pstd}Predicted probabilities of discrete choices for the median values of independent variables{p_end}
        . ziopprobabilities pr_choice

{pstd}Predicted probabilities of discrete choices for the specified values of independent variables{p_end}
        . ziopprobabilities, at (pb=1, spread=0.426, houst=1.6, gdp=6.8)

{pstd}Predicted probabilities of three types of zeros conditional on the regime for the median values of independent variables{p_end}
        . ziopprobabilities pr_zero, zeros

{pstd}Predicted probabilities of three regimes for the median values of independent variables{p_end}
        . ziopprobabilities pr_regime, regimes


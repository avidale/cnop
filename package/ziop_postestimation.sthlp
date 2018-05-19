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
{synopt :{cmd:zeros}}indicates that the probabilities of different types of zeros (the outcomes in the inflated category), conditional on different regimes, must be predicted instead of the choice probabilities{p_end}
{synopt :{cmd:regimes}}indicates that the probabilities of the regimes must be predicted instead of the choice probabilities; this option is ignored if the option {cmd:zeros}} is used.{p_end}
{synopt :{opt output(string)}}specifies the different types of predictions. The possible options for {it:string} are: {it:choice} for reporting the predicted outcome (the choice with the largest predicted probability); {it:mean} for reporting the expected value of the dependent variable computed as a summation of i*Pr(y=i) across all choices i; and 
{it:cum} for predicting the cumulative choice probabilities such as Pr(y<=0), Pr(y<=1), ... . If {it:string} is not specified, the usual choice probabilities such as Pr(y=0), Pr(y=1), ... are predicted and saved into new variables with the {varname} prefix.{p_end}
{synoptline}
{p2colreset}{...}

{title:Description for predict}

{pstd}
    {cmd:predict} creates new variables containing predictions such as the predicted probabilities of discrete choices, the predicted probabilities of the regimes or the types of zeros conditional on the regime, the expected values of the dependent variable for all sample values of the independent variables.


{marker ziopmargins}{...}
{title:Syntax for ziopmargins}

TBD

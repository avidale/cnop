{smcl}
{* *! version 0.0.1  08may2018}{...}
{title:Title}

Postestimation tools for {help ziop:zero-inflated ordered probit models} {cmd:nop}, {cmd:ziop2} and {cmd:ziop3}

{title:Description}

{pstd}
The following postestimation commands are available after {cmd:nop}, {cmd:ziop2} and {cmd:ziop3}: 

{synoptset 20 notes}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
{synopt :{helpb ziop postestimation##predict:predict}}predictions, residuals, influence statistics, and other diagnostic measures{p_end}
{synopt :{helpb ziop postestimation##ziopprobabilities:ziopprobabilities}}predicted probabilities at a single point{p_end}
{synopt :{helpb ziop postestimation##ziopcontrasts:ziopcontrasts}}differences in predicted probabilities{p_end}
{synopt :{helpb ziop postestimation##ziopmargins:ziopmargins}}marginal effects{p_end}
{synopt :{helpb ziop postestimation##ziopclassification:ziopclassification}}classification table and goodness-of-fit statistics{p_end}
{synopt :{helpb ziop postestimation##ziopvuong:ziopvuong}}non-nested Vuong test{p_end}
{synoptline}
{p2colreset}{...}


{marker predict}{...}
{title:Syntax for predict}

The predict command after the {cmd:nop}, {cmd:ziop2} and {cmd:ziop3} estimation commands produces either
predicted probabilities or "central" values of the responses.

{cmd:predict} {opth name(varname)} {ifin} [, {opt zeros} {opt regime} {opt output(string)} {opt at(string)}]

{opt name} is the name of predicted variable, if it is single, 
or prefix for names, if there are several predicted variables

{opt zeros} indicates that different types of zeros (i.e. "intrinsic zeros", or "positive zeros", or
"negative zeros") must be predicted instead of different response values.

{opt regime} indicates that different groups of response (negative, positive or zero) must be
predicted instead of different response values. This option is ignored if {opt zeros} option is on.

{opt output} specifies type of aggregating predicted probabilities of different response. 
Possible values are: 
    {opt mode} for reporting the outcome with the highest predicted probability, 
    {opt mean} for predicting the expected outcome, 
    {opt cum}  for predicting cumulative response probabilities. 
If not specified, raw response probabilities are predicted and placed into multiple variables with prefix {opt name}.


{marker ziopmargins}{...}
{title:Syntax for ziopmargins}

TBD

* Stata .do file to generate data with an increasing trend for the exposure of
* interest. Used to test the null-effect SR estimation called by the
* sr_wrapper() function

clear
set obs 1000000

local expo_rate = 1/(50 * 365)
local outc_rate = 1/(1000 * 365)
local max_fu = 365 * 5

survsim t_outc fail, ///
	distribution(gompertz) lambda(`outc_rate') gamma(0.003) maxtime(`max_fu')
replace t_outc = . if !fail

survsim t_expo expo, ///
	distribution(exponential) lambda(`expo_rate') maxtime(`max_fu')
replace t_expo = . if !expo

export delimited input/trend-data.csv, replace
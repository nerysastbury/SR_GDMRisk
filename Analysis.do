**************Finding change in intervention and control groups **************
clear all
use "Redacted\trimmeddataset.dta"

************************************Clustering***********************************************************************************************************************
//Ferrara et al is a cluster randomised study - I need to account for this clustering my estimation of the SD. I need to do this for both baseline and follow-up N. 
	// they reported using a range of ICC values of 0.01-0.05? I will take the halfway point - 0.025, as the ICC here. 
	// to calculate the effective sample size as shown in the Cochrane handbook, I need to first work out the average cluster size 
	// this is (intv N + ctrl N)/(intv clusters + ctrl clusters) = (1087 + 1193)/(22+22) = 51.8 = 51.8. 
	// Now we need to calculate the design effect - this is 1+(average cluster size -1)*ICC = 1+(51.8-1)*0.025 = 2.27
	// so the INTV effective sample size = 1087/2.27 = 478.85 = 479, and for the CTRL, it's  1193/2.27 = 525.55 = 526
replace intv_baseline_n = 479 if study=="Ferrara, 2016"
replace ctrl_baseline_n = 526 if study=="Ferrara, 2016"

	//now we need to do this for follow-up N
	// average no. of patients per cluster: intv = 676, ctrl = 744 => (676 + 744)/(44) = 32.3
	//design effect = 1+(32.3-1)*0.025 = 1.78
	//so intv follow up N = 676/1.78 = 379.8 = 380, and ctrl follow up N = 417.9 = 418
	
replace intv_fup_n = 380 if study=="Ferrara, 2016"
replace ctrl_fup_n = 418 if study=="Ferrara, 2016"

	//and change n in the table to the effective sample size 
replace n =1005 if study=="Ferrara, 2016"

//as for T2DM count for Ferrara, Cochrane says: "For dichotomous data, both the number of participants and the number experiencing the event should be divided by the same design effect."
	 //so:
	 replace intv_fup_mean = intv_fup_mean/1.78 if study=="Ferrara, 2016" & bio=="t2dm"
	 replace ctrl_fup_mean = ctrl_fup_mean/1.78 if study=="Ferrara, 2016" & bio=="t2dm"

*******************Let's do some dichotomous data analysis**********************

////create the 2x2 cells: a = events intervention, b = non-events intervention,
*                          c = events control,     d = non-events control
generate a = intv_fup_mean if bio=="t2dm"
generate b = intv_fup_n - intv_fup_mean if bio=="t2dm"
generate c = ctrl_fup_mean if bio=="t2dm"
generate d = ctrl_fup_n - ctrl_fup_mean if bio=="t2dm"

* keep backups (so you can revert if you try CC later)
clonevar a_orig = a
clonevar b_orig = b
clonevar c_orig = c
clonevar d_orig = d

* --- 3. flag zero-event studies for inspection
generate any_zero_event = (a==0 | c==0)
tab any_zero_event
list study a b c d any_zero_event if any_zero_event==1


* --- 4. compute RR, logRR and SE(logRR)
* Note: logRR and se will be missing / infinite if any event cell is zero.
generate rr = . 
replace rr = (a/(a+b)) / (c/(c+d)) if (a>=0 & c>=0 & (a+b)>0 & (c+d)>0)

generate logrr = .
replace logrr = ln(rr) if rr>0

generate se_logrr = .
replace se_logrr = sqrt( (1/a) - (1/(a+b)) + (1/c) - (1/(c+d)) ) ///
    if a>0 & c>0
	
save "Redacted\trimmeddataset.dta", replace


******************************************************************************************************************
****************************** r = 0.8*******************************************
******************************************************************************************************************

**************Finding change in intervention and control groups **************
clear all
use "Redacted\trimmeddataset.dta"

//////////////////// find difference for intervention
gen intv_change_mean=.
replace intv_change_mean = intv_fup_mean if changevsactual=="change"
replace intv_change_mean = intv_fup_mean - intv_base_mean if changevsactual=="actual"

gen intv_change_sd=.
replace intv_change_sd = intv_fup_var if changevsactual=="change"

//HERE IS THE RIGHT WAY - I will start by assuming an r of 0.8, and then I will use other r values for sensitivity analyses - so I will create a separate dataset for each of these 
local r = 0.8
replace intv_change_sd = sqrt( intv_base_var^2 + intv_fup_var^2 -( 2 * `r' * intv_base_var * intv_fup_var )) if changevsactual == "actual"

//////////////////////////// find difference for control
gen ctrl_change_mean=.
replace ctrl_change_mean = ctrl_fup_mean if controlactualvschance=="change"
replace ctrl_change_mean = ctrl_fup_mean - ctrl_base_mean if controlactualvschance=="actual"

gen ctrl_change_sd=.
replace ctrl_change_sd = ctrl_fup_var if controlactualvschance=="change"
//I am going to set my correlation coefficient at 0.8
local r = 0.8
replace ctrl_change_sd = sqrt( ctrl_base_var^2 + ctrl_fup_var^2 -( 2 * `r' * ctrl_base_var * ctrl_fup_var )) if controlactualvschance == "actual"


save "Redacted\corr0.8.dta", replace

******************************************************************BOCFBOCFBOCF*******************************************************************
**************and now we do baseline carried forward (bocf) analysis for intervention grp
clear all
use "Redacted\corr0.8.dta"

gen intv_change_mean_obs = intv_change_mean
gen ctrl_change_mean_obs = ctrl_change_mean

replace intv_change_mean = (intv_change_mean*intv_fup_n)/intv_baseline_n
//now get variance
replace  intv_change_sd = (intv_change_sd^2*(intv_fup_n -1) + ((intv_fup_n/intv_baseline_n)*(1-(intv_fup_n/intv_baseline_n))*intv_baseline_n*intv_change_mean_obs^2))/(intv_baseline_n -1)
//now get sd 
replace intv_change_sd = sqrt(intv_change_sd)

gen intv_change_se=.
replace intv_change_se = intv_change_sd/sqrt(intv_baseline_n)

**************and now we do baseline observation carried forward (bocf) analysis for control grp
replace ctrl_change_mean = (ctrl_change_mean*ctrl_fup_n)/ctrl_baseline_n
//now get variance
replace  ctrl_change_sd = (ctrl_change_sd^2*(ctrl_fup_n -1) + ((ctrl_fup_n/ctrl_baseline_n)*(1-(ctrl_fup_n/ctrl_baseline_n))*ctrl_baseline_n*ctrl_change_mean_obs^2))/(ctrl_baseline_n -1)
//now get sd 
replace ctrl_change_sd = sqrt(ctrl_change_sd)

gen ctrl_change_se=.
replace ctrl_change_se = ctrl_change_sd/sqrt(ctrl_baseline_n)

*********************************Finding mean difference********************
gen mean_difference=.
replace mean_difference = intv_change_mean - ctrl_change_mean
 
gen mean_difference_se=.
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2))
	

//but then for hba1c mcmanus and some casey gives post-intervention scores rather than change from baseline, so I will (1) change sds into ses, (2) do intv mean - ctrl mean, then (3) use the same se ///
///formula as the above 

//so first let's convert sds to ses 
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="bmi"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="wc"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"


replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="bmi"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="wc"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"

//now calculate change se
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="bmi"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="wc"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"

//now do intv-ctrl 
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="McManus, 2018" & bio=="HbA1c (%)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="bmi"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="wc"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"


//and for quansah and Highton I need to input the mean difference and se I calculated from the cis they reported, because they don't report actual numbers at baseline-follow-up. I calculated SE from 9%% CIs as (upper interval - lower interval/3.92)
replace mean_difference = -0.38 if study=="Quansah, 2024" & bio=="weight"
replace mean_difference = -0.09 if study=="Quansah, 2024" & bio=="Fasting plasma glucose (mmol/l)"
replace mean_difference = 0.21858 if study=="Quansah, 2024" & bio=="HbA1c (mmol/mol)"
replace mean_difference = 0.02 if study=="Quansah, 2024" & bio=="HbA1c (%)"
replace mean_difference = -0.09 if study=="Quansah, 2024" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference = -0.02 if study=="Khunti, 2022" & bio=="High-density lipoprotein cholesterol (mmol/l)"
replace mean_difference = 0.01 if study=="Khunti, 2022" & bio=="Low-density lipoprotein cholesterol (mmol/l)"
replace mean_difference = -0.03 if study=="Khunti, 2022" & bio=="Total cholesterol (mmol/l)"
replace mean_difference = 0.09 if study=="Khunti, 2022" & bio=="Triglycerides (mmol/l)"
replace mean_difference = 0.89 if study=="Khunti, 2022" & bio=="Systolic blood pressure (mmHg)"
replace mean_difference = -0.4 if study=="Khunti, 2022" & bio=="weight"
replace mean_difference = -0.48 if study=="Khunti, 2022" & bio=="wc"
replace mean_difference = -0.17 if study=="Khunti, 2022" & bio=="HbA1c (mmol/mol)"
replace mean_difference = -0.04 if study=="Khunti, 2022" & bio=="HbA1c (%)"
replace mean_difference = -0.18 if study=="Khunti, 2022" & bio=="bmi"
replace mean_difference = -0.56 if study=="Khunti, 2022" & bio=="Diastolic blood pressure (mmHg)"

//the SEs below were calculated using SE = (upper CI - lower CI)/(3.92) for 95% CIs
replace mean_difference_se = ((1.30 - -2.08)/3.92) if study=="Quansah, 2024" & bio=="weight"
replace mean_difference_se = ((0.08 - -0.26)/3.92) if study=="Quansah, 2024" & bio=="Fasting plasma glucose (mmol/l)"
replace mean_difference_se = ((((0.09*10.929)-23.5)- ((-0.05*10.929)-23.5))/3.92) if study=="Quansah, 2024" & bio=="HbA1c (mmol/mol)" //I have to convert CIs to mmol/mol and calculate SE)
replace mean_difference_se = ((0.09 - -0.05)/3.92) if study=="Quansah, 2024" & bio=="HbA1c (%)"
replace mean_difference_se = ((0.55 - -0.53)/3.92) if study=="Quansah, 2024" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference_se = ((0.06 - -0.09)/3.92) if study=="Khunti, 2022" & bio=="High-density lipoprotein cholesterol (mmol/l)"
replace mean_difference_se = ((0.14 - -0.12)/3.92) if study=="Khunti, 2022" & bio=="Low-density lipoprotein cholesterol (mmol/l)"
replace mean_difference_se = ((0.14 - -0.19)/3.92) if study=="Khunti, 2022" & bio=="Total cholesterol (mmol/l)"
replace mean_difference_se = ((0.28 - -0.09)/3.92) if study=="Khunti, 2022" & bio=="Triglycerides (mmol/l)"
replace mean_difference_se = ((3.30 - -1.52)/3.92) if study=="Khunti, 2022" & bio=="Systolic blood pressure (mmHg)"
replace mean_difference_se = ((1.04 - -1.83)/3.92) if study=="Khunti, 2022" & bio=="weight"
replace mean_difference_se = ((1.63 - -2.58)/3.92) if study=="Khunti, 2022" & bio=="wc"
replace mean_difference_se = ((0.68 - -1.03)/3.92) if study=="Khunti, 2022" & bio=="HbA1c (mmol/mol)"
replace mean_difference_se = ((0.04 - -0.012)/3.92) if study=="Khunti, 2022" & bio=="HbA1c (%)"
replace mean_difference_se = ((0.36 - -0.73 )/3.92) if study=="Khunti, 2022" & bio=="bmi"
replace mean_difference_se = ((1.25 - -2.37)/3.92) if study=="Khunti, 2022" & bio=="Diastolic blood pressure (mmHg)"
save "Redacted\BCOFr0.8.dta", replace


//export to excel
cd "Redacted"
export excel using "BCOFanalysisr0.8xls", firstrow(variables) replace
//
//
//
//
//
//
//
******************************************************************************************************************
****************************** r = 0.6*******************************************
******************************************************************************************************************
//
//
//
//
//
//
//
//
//
**************Finding change in intervention and control groups **************
clear all
use "Redacted\trimmeddataset.dta"

//////////////////// find difference for intervention
gen intv_change_mean=.
replace intv_change_mean = intv_fup_mean if changevsactual=="change"
replace intv_change_mean = intv_fup_mean - intv_base_mean if changevsactual=="actual"

gen intv_change_sd=.
replace intv_change_sd = intv_fup_var if changevsactual=="change"


//HERE IS THE RIGHT WAY - I will start by assuming an r of 0.6, and then I will use other r values for sensitivity analyses - so I will create a separate dataset for each of these 
local r = 0.6
replace intv_change_sd = sqrt( intv_base_var^2 + intv_fup_var^2 -( 2 * `r' * intv_base_var * intv_fup_var )) if changevsactual == "actual"

//////////////////////////// find difference for control
gen ctrl_change_mean=.
replace ctrl_change_mean = ctrl_fup_mean if controlactualvschance=="change"
replace ctrl_change_mean = ctrl_fup_mean - ctrl_base_mean if controlactualvschance=="actual"

gen ctrl_change_sd=.
replace ctrl_change_sd = ctrl_fup_var if controlactualvschance=="change"
//I am going to set my correlation coefficient at 0.6
local r = 0.6
replace ctrl_change_sd = sqrt( ctrl_base_var^2 + ctrl_fup_var^2 -( 2 * `r' * ctrl_base_var * ctrl_fup_var )) if controlactualvschance == "actual"


save "Redacted\corr0.6.dta", replace

******************************************************************BOCFBOCFBOCF*******************************************************************
**************and now we do baseline carried forward (bocf) analysis for intervention grp
clear all
use "Redacted\corr0.6.dta"

gen intv_change_mean_obs = intv_change_mean
gen ctrl_change_mean_obs = ctrl_change_mean

replace intv_change_mean = (intv_change_mean*intv_fup_n)/intv_baseline_n
//now get variance
replace  intv_change_sd = (intv_change_sd^2*(intv_fup_n -1) + ((intv_fup_n/intv_baseline_n)*(1-(intv_fup_n/intv_baseline_n))*intv_baseline_n*intv_change_mean_obs^2))/(intv_baseline_n -1)
//now get sd 
replace intv_change_sd = sqrt(intv_change_sd)

gen intv_change_se=.
replace intv_change_se = intv_change_sd/sqrt(intv_baseline_n)

**************and now we do baseline observation carried forward (bocf) analysis for control grp
replace ctrl_change_mean = (ctrl_change_mean*ctrl_fup_n)/ctrl_baseline_n
//now get variance
replace  ctrl_change_sd = (ctrl_change_sd^2*(ctrl_fup_n -1) + ((ctrl_fup_n/ctrl_baseline_n)*(1-(ctrl_fup_n/ctrl_baseline_n))*ctrl_baseline_n*ctrl_change_mean_obs^2))/(ctrl_baseline_n -1)
//now get sd 
replace ctrl_change_sd = sqrt(ctrl_change_sd)

gen ctrl_change_se=.
replace ctrl_change_se = ctrl_change_sd/sqrt(ctrl_baseline_n)

*********************************Finding mean difference********************
gen mean_difference=.
replace mean_difference = intv_change_mean - ctrl_change_mean
 
gen mean_difference_se=.
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2))
	

//but then for hba1c mcmanus and some casey gives post-intervention scores rather than change from baseline, so I will (1) change sds into ses, (2) do intv mean - ctrl mean, then (3) use the same se ///
///formula as the above 

//so first let's convert sds to ses 
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="bmi"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="wc"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"


replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="bmi"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="wc"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"

//now calculate change se
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="bmi"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="wc"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"

//now do intv-ctrl 
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="McManus, 2018" & bio=="HbA1c (%)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="bmi"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="wc"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"


//and for quansah and Highton I need to input the mean difference and se I calculated from the cis they reported, because they don't report actual numbers at baseline-follow-up. I calculated SE from 9%% CIs as (upper interval - lower interval/3.92)
replace mean_difference = -0.38 if study=="Quansah, 2024" & bio=="weight"
replace mean_difference = -0.09 if study=="Quansah, 2024" & bio=="Fasting plasma glucose (mmol/l)"
replace mean_difference = 0.21858 if study=="Quansah, 2024" & bio=="HbA1c (mmol/mol)"
replace mean_difference = 0.02 if study=="Quansah, 2024" & bio=="HbA1c (%)"
replace mean_difference = -0.09 if study=="Quansah, 2024" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference = -0.02 if study=="Khunti, 2022" & bio=="High-density lipoprotein cholesterol (mmol/l)"
replace mean_difference = 0.01 if study=="Khunti, 2022" & bio=="Low-density lipoprotein cholesterol (mmol/l)"
replace mean_difference = -0.03 if study=="Khunti, 2022" & bio=="Total cholesterol (mmol/l)"
replace mean_difference = 0.09 if study=="Khunti, 2022" & bio=="Triglycerides (mmol/l)"
replace mean_difference = 0.89 if study=="Khunti, 2022" & bio=="Systolic blood pressure (mmHg)"
replace mean_difference = -0.4 if study=="Khunti, 2022" & bio=="weight"
replace mean_difference = -0.48 if study=="Khunti, 2022" & bio=="wc"
replace mean_difference = -0.17 if study=="Khunti, 2022" & bio=="HbA1c (mmol/mol)"
replace mean_difference = -0.04 if study=="Khunti, 2022" & bio=="HbA1c (%)"
replace mean_difference = -0.18 if study=="Khunti, 2022" & bio=="bmi"
replace mean_difference = -0.56 if study=="Khunti, 2022" & bio=="Diastolic blood pressure (mmHg)"

//the SEs below were calculated using SE = (upper CI - lower CI)/(3.92) for 95% CIs
replace mean_difference_se = ((1.30 - -2.08)/3.92) if study=="Quansah, 2024" & bio=="weight"
replace mean_difference_se = ((0.08 - -0.26)/3.92) if study=="Quansah, 2024" & bio=="Fasting plasma glucose (mmol/l)"
replace mean_difference_se = ((((0.09*10.929)-23.5)- ((-0.05*10.929)-23.5))/3.92) if study=="Quansah, 2024" & bio=="HbA1c (mmol/mol)" //I have to convert CIs to mmol/mol and calculate SE)
replace mean_difference_se = ((0.09 - -0.05)/3.92) if study=="Quansah, 2024" & bio=="HbA1c (%)"
replace mean_difference_se = ((0.55 - -0.53)/3.92) if study=="Quansah, 2024" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference_se = ((0.06 - -0.09)/3.92) if study=="Khunti, 2022" & bio=="High-density lipoprotein cholesterol (mmol/l)"
replace mean_difference_se = ((0.14 - -0.12)/3.92) if study=="Khunti, 2022" & bio=="Low-density lipoprotein cholesterol (mmol/l)"
replace mean_difference_se = ((0.14 - -0.19)/3.92) if study=="Khunti, 2022" & bio=="Total cholesterol (mmol/l)"
replace mean_difference_se = ((0.28 - -0.09)/3.92) if study=="Khunti, 2022" & bio=="Triglycerides (mmol/l)"
replace mean_difference_se = ((3.30 - -1.52)/3.92) if study=="Khunti, 2022" & bio=="Systolic blood pressure (mmHg)"
replace mean_difference_se = ((1.04 - -1.83)/3.92) if study=="Khunti, 2022" & bio=="weight"
replace mean_difference_se = ((1.63 - -2.58)/3.92) if study=="Khunti, 2022" & bio=="wc"
replace mean_difference_se = ((0.68 - -1.03)/3.92) if study=="Khunti, 2022" & bio=="HbA1c (mmol/mol)"
replace mean_difference_se = ((0.04 - -0.012)/3.92) if study=="Khunti, 2022" & bio=="HbA1c (%)"
replace mean_difference_se = ((0.36 - -0.73 )/3.92) if study=="Khunti, 2022" & bio=="bmi"
replace mean_difference_se = ((1.25 - -2.37)/3.92) if study=="Khunti, 2022" & bio=="Diastolic blood pressure (mmHg)"
save "Redacted\BCOFr0.6.dta", replace


//export to excel
cd "Redacted"
export excel using "BCOFanalysisr0.6xls", firstrow(variables) replace


//
//
//
//
//
//
//
//
//
******************************************************************************************************************
****************************** r = 0.8 without BOCF*******************************************
******************************************************************************************************************
//
//
//
//
//
//
//
**************Finding change in intervention and control groups **************
clear all
use "Redacted\trimmeddataset.dta"

//////////////////// find difference for intervention
gen intv_change_mean=.
replace intv_change_mean = intv_fup_mean if changevsactual=="change"
replace intv_change_mean = intv_fup_mean - intv_base_mean if changevsactual=="actual"

gen intv_change_sd=.
replace intv_change_sd = intv_fup_var if changevsactual=="change"

//BELOW IS THE OLD *WRONG* WAY TO CALCULATE CHANGE SD
//replace intv_change_sd = 0.96 * sqrt(intv_base_var^2) * sqrt(intv_fup_var^2) if changevsactual=="actual"
//replace intv_change_sd = sqrt(intv_base_var^2 + intv_fup_var^2 -(2*intv_change_sd)) if changevsactual=="actual"

//HERE IS THE RIGHT WAY - I will start by assuming an r of 0.8, and then I will use other r values for sensitivity analyses - so I will create a separate dataset for each of these 
local r = 0.8
replace intv_change_sd = sqrt( intv_base_var^2 + intv_fup_var^2 -( 2 * `r' * intv_base_var * intv_fup_var )) if changevsactual == "actual"

//////////////////////////// find difference for control
gen ctrl_change_mean=.
replace ctrl_change_mean = ctrl_fup_mean if controlactualvschance=="change"
replace ctrl_change_mean = ctrl_fup_mean - ctrl_base_mean if controlactualvschance=="actual"

gen ctrl_change_sd=.
replace ctrl_change_sd = ctrl_fup_var if controlactualvschance=="change"
//I am going to set my correlation coefficient at 0.8
local r = 0.8
replace ctrl_change_sd = sqrt( ctrl_base_var^2 + ctrl_fup_var^2 -( 2 * `r' * ctrl_base_var * ctrl_fup_var )) if controlactualvschance == "actual"

gen intv_change_se=.
replace intv_change_se = intv_change_sd/sqrt(intv_baseline_n)

gen ctrl_change_se=.
replace ctrl_change_se = ctrl_change_sd/sqrt(ctrl_baseline_n)

*********************************Finding mean difference********************
gen mean_difference=.
replace mean_difference = intv_change_mean - ctrl_change_mean
 
gen mean_difference_se=.
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2))
	

//but then for hba1c mcmanus and some casey gives post-intervention scores rather than change from baseline, so I will (1) change sds into ses, (2) do intv mean - ctrl mean, then (3) use the same se ///
///formula as the above 

//so first let's convert sds to ses 
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="bmi"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="wc"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace intv_change_se = intv_fup_var/sqrt(intv_fup_n) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"


replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="bmi"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="wc"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace ctrl_change_se = ctrl_fup_var/sqrt(ctrl_fup_n) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"

//now calculate change se
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="McManus, 2018" & bio=="HbA1c (%)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="bmi"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="wc"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference_se = sqrt((intv_change_se^2)+(ctrl_change_se^2)) if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"

//now do intv-ctrl 
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="McManus, 2018" & bio=="HbA1c (mmol/mol)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="McManus, 2018" & bio=="HbA1c (%)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="bmi"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="wc"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference = intv_fup_mean - ctrl_fup_mean if study=="Casey, 2020" & bio=="Fasting plasma glucose (mmol/l)"


//and for quansah and Highton I need to input the mean difference and se I calculated from the cis they reported, because they don't report actual numbers at baseline-follow-up. I calculated SE from 9%% CIs as (upper interval - lower interval/3.92)
replace mean_difference = -0.38 if study=="Quansah, 2024" & bio=="weight"
replace mean_difference = -0.09 if study=="Quansah, 2024" & bio=="Fasting plasma glucose (mmol/l)"
replace mean_difference = 0.21858 if study=="Quansah, 2024" & bio=="HbA1c (mmol/mol)"
replace mean_difference = 0.02 if study=="Quansah, 2024" & bio=="HbA1c (%)"
replace mean_difference = -0.09 if study=="Quansah, 2024" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference = -0.02 if study=="Khunti, 2022" & bio=="High-density lipoprotein cholesterol (mmol/l)"
replace mean_difference = 0.01 if study=="Khunti, 2022" & bio=="Low-density lipoprotein cholesterol (mmol/l)"
replace mean_difference = -0.03 if study=="Khunti, 2022" & bio=="Total cholesterol (mmol/l)"
replace mean_difference = 0.09 if study=="Khunti, 2022" & bio=="Triglycerides (mmol/l)"
replace mean_difference = 0.89 if study=="Khunti, 2022" & bio=="Systolic blood pressure (mmHg)"
replace mean_difference = -0.4 if study=="Khunti, 2022" & bio=="weight"
replace mean_difference = -0.48 if study=="Khunti, 2022" & bio=="wc"
replace mean_difference = -0.17 if study=="Khunti, 2022" & bio=="HbA1c (mmol/mol)"
replace mean_difference = -0.04 if study=="Khunti, 2022" & bio=="HbA1c (%)"
replace mean_difference = -0.18 if study=="Khunti, 2022" & bio=="bmi"
replace mean_difference = -0.56 if study=="Khunti, 2022" & bio=="Diastolic blood pressure (mmHg)"

//the SEs below were calculated using SE = (upper CI - lower CI)/(3.92) for 95% CIs
replace mean_difference_se = ((1.30 - -2.08)/3.92) if study=="Quansah, 2024" & bio=="weight"
replace mean_difference_se = ((0.08 - -0.26)/3.92) if study=="Quansah, 2024" & bio=="Fasting plasma glucose (mmol/l)"
replace mean_difference_se = ((((0.09*10.929)-23.5)- ((-0.05*10.929)-23.5))/3.92) if study=="Quansah, 2024" & bio=="HbA1c (mmol/mol)" //I have to convert CIs to mmol/mol and calculate SE)
replace mean_difference_se = ((0.09 - -0.05)/3.92) if study=="Quansah, 2024" & bio=="HbA1c (%)"
replace mean_difference_se = ((0.55 - -0.53)/3.92) if study=="Quansah, 2024" & bio=="2-hour Oral glucose tolerance test (mmol/l)"
replace mean_difference_se = ((0.06 - -0.09)/3.92) if study=="Khunti, 2022" & bio=="High-density lipoprotein cholesterol (mmol/l)"
replace mean_difference_se = ((0.14 - -0.12)/3.92) if study=="Khunti, 2022" & bio=="Low-density lipoprotein cholesterol (mmol/l)"
replace mean_difference_se = ((0.14 - -0.19)/3.92) if study=="Khunti, 2022" & bio=="Total cholesterol (mmol/l)"
replace mean_difference_se = ((0.28 - -0.09)/3.92) if study=="Khunti, 2022" & bio=="Triglycerides (mmol/l)"
replace mean_difference_se = ((3.30 - -1.52)/3.92) if study=="Khunti, 2022" & bio=="Systolic blood pressure (mmHg)"
replace mean_difference_se = ((1.04 - -1.83)/3.92) if study=="Khunti, 2022" & bio=="weight"
replace mean_difference_se = ((1.63 - -2.58)/3.92) if study=="Khunti, 2022" & bio=="wc"
replace mean_difference_se = ((0.68 - -1.03)/3.92) if study=="Khunti, 2022" & bio=="HbA1c (mmol/mol)"
replace mean_difference_se = ((0.04 - -0.012)/3.92) if study=="Khunti, 2022" & bio=="HbA1c (%)"
replace mean_difference_se = ((0.36 - -0.73 )/3.92) if study=="Khunti, 2022" & bio=="bmi"
replace mean_difference_se = ((1.25 - -2.37)/3.92) if study=="Khunti, 2022" & bio=="Diastolic blood pressure (mmHg)"
save "Redacted\noBCOFr0.8.dta", replace


//export to excel
cd "Redacted"
export excel using "noBCOFanalysisr0.8xls", firstrow(variables) replace


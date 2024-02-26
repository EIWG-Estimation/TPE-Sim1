/*******************************************************************
| Name       : analysis_plots.sas
| Purpose    : Plots from the summary datasets 
|              and macros. 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 15SEP21 
|-------------------------------------------------------------------
| Notes:
|
| 
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";

proc format lib = work;
  value rdrate 1 = "1"
               2 = "2"
               3 = "3"
               4 = "4"
               5 = "5"
               6 = "6";
       
  value scenario 1 = "DAR Instant Change"
                 2 = "DAR Gradual Change"
                 3 = "DNAR Instant Change"
                 4 = "DNAR Gradual Change";                    
run;


********************************************************************;
*** STACK SUMMARY DATA                                           ***;
********************************************************************;

data all_dif;
  set results.dar_drop_summary_dif_all       (in = in1)
      results.dar_linearchg_summary_dif_all  (in = in2)
      results.dnar_drop_summary_dif_all      (in = in3)
      results.dnar_linearchg_summary_dif_all (in = in4);
;
  if in1 then scenario = 1;
  else if in2 then scenario = 2;
  else if in3 then scenario = 3;
  else if in4 then scenario = 4;
   
  dif_policy_se_mean1 = dif_policy_se_mean;
  if method in ("MMRM2" "MMRM3") then dif_policy_se_mean2 = dif_policy_se_adj_mean;
  else dif_policy_se_mean2 = dif_policy_se_mean;  

  dif_policy_rmse_mean1 = dif_policy_rmse_mean;
  if method in ("MMRM2" "MMRM3") then dif_policy_rmse_mean2 = dif_policy_rmse_adj_mean;
  else dif_policy_rmse_mean2 = dif_policy_rmse_mean;  
   
  dif_policy_hw_mean1 = dif_policy_hw_mean;
  if method in ("MMRM2" "MMRM3") then dif_policy_hw_mean2 = dif_policy_hw_adj_mean;
  else dif_policy_hw_mean2 = dif_policy_hw_mean;

  dif_policy_cic_mean1 = dif_policy_cic_mean;
  if method in ("MMRM2" "MMRM3") then dif_policy_cic_mean2 = dif_policy_cic_adj_mean;
  else dif_policy_cic_mean2 = dif_policy_cic_mean;
  
run;


********************************************************************;
*** CREATE ATTRIBUTE MAP FOR METHODS  #8F00FF                        ***;
********************************************************************;
  
data map1;
  array values   [10] $10 ("FULL"  "MMRM1"  "MMRM2"  "MMRM3"  "MI1"  "MI2"  "MI3" "JTR" "CIR" "CR");
  array colors   [10] $10 ("#000000"  "#88C6ED" "Pink" "#90EE90"  "#394BA0" "Darkred" "DarkGreen" "#54d2d2" "#ffcb00" "#ff6150");
  array patterns [10] $10 ("34"  "1"  "1"  "1"  "20"  "20"  "20"  "1"  "1"  "1");
  id = "models";
  do i = 1 to dim(values);
    value       = values[i];
    linecolor   = colors[i];
    linepattern = patterns[i];  
    linethickness = 2;
    output;
  end;
run;  


********************************************************************;
*** LINE PLOTS FOR BIAS "#CC99FF"                                ***;
********************************************************************;

ods graphics on / width = 8in height=6in; 

*** RD METHODS ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 1 2 3 4 5 6) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_bias_mean / group = method attrid=models;
  title1 "Figure 1A - Treatment effect bias";
  rowaxis label = "Bias" values=(-0.10 to 0.04 by 0.02) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;

*** RBI METHODS ***;  
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 7 8 9))) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_bias_mean / group = method attrid=models;
  title1 "Figure 1B - Treatment effect bias";
  rowaxis label = "Bias" values=(-0.10 to 0.04 by 0.02) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;


********************************************************************;
*** LINE PLOTS FOR SD(POINT ESTIMATE)                            ***;
********************************************************************;

*** RD METHODS ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 1 2 3 4 5 6) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_est_sd / group = method attrid=models;
  title1 "Figure 2A - SD of treatement effect estimates";
  rowaxis label = "SD of treatment effect estimate" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;

*** RBI METHODS ***;  
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 7 8 9))) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_est_sd / group = method attrid=models;
  title1 "Figure 2B - SD of treatement effect estimates";
  rowaxis label = "SD of Treatment Effect Estimate" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;



********************************************************************;
*** LINE PLOTS FOR AVERAGE MODEL SE                              ***;
********************************************************************;

*** RD METHODS - UNADJUSTED ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 1 2 3 4 5 6) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_se_mean1 / group = method attrid=models;
  title1 "Figure 3A - Treatment Effect Model SE (unadjusted*)";
  footnote1 "* No variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Average Model SE" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;

*** RD METHODS - ADJUSTED ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 1 2 3 4 5 6) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_se_mean2 / group = method attrid=models;
  title1 "Figure 3B - Treatment Effect Model SE (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Average Model SE" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;

*** RBI METHODS ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 7 8 9) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_se_mean1 / group = method attrid=models;
  title1 "Figure 3C - Treatement effect Model SE";
  rowaxis label = "Average Model SE" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;


********************************************************************;
*** LINE PLOTS FOR ROOT OF AVERAGE MSE                           ***;
********************************************************************;

*** RD METHODS - UNADJUSTED ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 1 2 3 4 5 6) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_rmse_mean1 / group = method attrid=models;
  title1 "Figure 5A - Treatment Effect Root of Average MSE (unadjusted*)";
  footnote1 "* No variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Root of Average MSE" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;

*** RD METHODS - ADJUSTED ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 1 2 3 4 5 6) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_se_mean2 / group = method attrid=models;
  title1 "Figure 5B - Treatment Effect Root of Average MSE (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Root of Average MSE" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;

*** RBI METHODS ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 7 8 9) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_se_mean1 / group = method attrid=models;
  title1 "Figure 5C - Treatment Effect Root of Average MSE";
  rowaxis label = "Root of Average MSE" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;


********************************************************************;
*** LINE PLOTS FOR HALFWIDTHS                                    ***;
********************************************************************;

*** RD METHODS - UNADJUSTED ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 1 2 3 4 5 6) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_hw_mean1 / group = method attrid=models;
  title1 "Treatement effect halfwidth (unadjusted*) vs Conditional missingness - RD methods";
  footnote1 "* No variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Halfwidth" values=(0.2 to 0.35 by 0.025) grid;
  colaxis label = "Conditional missingness rate" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;

*** RD METHODS - ADJUSTED ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 1 2 3 4 5 6) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_hw_mean2 / group = method attrid=models;
  title1 "Treatement effect halfwidth (adjusted*) vs Conditional missingness - RD methods";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Halfwidth" values=(0.2 to 0.35 by 0.025) grid;
  colaxis label = "Conditional missingness rate" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;

*** RBI METHODS ***;
proc sgpanel data = all_dif (where = (visitn = 5 and methodn in (0 7 8 9) )) dattrmap=map1;
  panelby scenario / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_hw_mean1 / group = method attrid=models;
  title1 "Treatement effect halfwidth vs Conditional missingness - RBI methods";
  rowaxis label = "Halfwidth" values=(0.2 to 0.35 by 0.025) grid;
  colaxis label = "Conditional missingness rate" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. scenario scenario.;
run;





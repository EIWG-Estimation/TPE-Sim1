/*******************************************************************
| Name       : analysis_plots_all.sas
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

*** SET UP HTML AND GRAPHICS FILE ROUTING ***;
filename out "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/output";

proc format lib = work;
  value rdrate 1 = "1"
               2 = "2"
               3 = "3"
               4 = "4"
               5 = "5"
               6 = "6";
       
  value typen 1 = "DAR Instant Change - Simple and RD Approaches"
              2 = "DAR Gradual Change - Simple and RD Approaches"
              3 = "DAR Instant Change - RBI Approaches"
              4 = "DAR Gradual Change - RBI Approaches"
              5 = "DNAR Instant Change - Simple and RD Approaches"
              6 = "DNAR Gradual Change - Simple and RD Approaches"
              7 = "DNAR Instant Change - RBI Approaches"
              8 = "DNAR Gradual Change - RBI Approaches"
             ;                    
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

  dif_policy_pwr_mean1 = dif_policy_pwr_mean;
  if method in ("MMRM2" "MMRM3") then dif_policy_pwr_mean2 = dif_policy_pwr_adj_mean;
  else dif_policy_pwr_mean2 = dif_policy_pwr_mean;

  dif_policy_pwr2_mean1 = dif_policy_pwr2_mean;
  if method in ("MMRM2" "MMRM3") then dif_policy_pwr2_mean2 = dif_policy_pwr2_adj_mean;
  else dif_policy_pwr2_mean2 = dif_policy_pwr2_mean;
 
  dif_policy_se_bias1 = dif_policy_se_mean1 - dif_policy_est_sd; 
  dif_policy_se_bias2 = dif_policy_se_mean2 - dif_policy_est_sd;
  
run;

data all_dif_plot;
  set all_dif (in = in1 where = (method in ("FULL"  "MMRM1"  "MMRM2"  "MMRM3"  "MI1"  "MI2"  "MI3")))
      all_dif (in = in2 where = (method in ("FULL"  "CIR" "CR" "JTR")));
  
  
  if in1 then do;
    if scenario = 1 then typen = 1;
    else if scenario = 2 then typen = 2;
    else if scenario = 3 then typen = 5;
    else if scenario = 4 then typen = 6;
  end;
  else if in2 then do;
    if scenario = 1 then typen = 3;
    else if scenario = 2 then typen = 4;
    else if scenario = 3 then typen = 7;
    else if scenario = 4 then typen = 8;
  end;
    
run;


********************************************************************;
*** CREATE ATTRIBUTE MAP FOR METHODS  #8F00FF                        ***;
********************************************************************;
  
data map1;
  array values   [10] $10 ("FULL"  "MMRM1"  "MMRM2"  "MMRM3"  "MI1"  "MI2"  "MI3" "CIR" "CR" "JTR");
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
*** MAIN FIGURES                                                 ***;
********************************************************************;

*** FIG 1 TREATMENT EFFECT BIAS - DAR SCENARIOS ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="FIG1";
ods html path=out gpath=out file="fig1.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_bias_mean / group = method attrid=models;
  title1 "Figure 1 - Treatment effect bias";
  rowaxis label = "Bias" values=(-0.10 to 0.04 by 0.02) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;


*** FIG 2 TREATMENT EFFECT SD(POINT ESTIMATE) ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="FIG2";
ods html path=out gpath=out file="fig2.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_est_sd / group = method attrid=models;
  title1 "Figure 2 - SD of treatement effect estimates";
  rowaxis label = "SD of treatment effect estimate" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;


*** FIG 3 TREATMENT EFFECT MODEL SE - ADJUSTED ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="FIG3";
ods html path=out gpath=out file="fig3.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_se_mean2 / group = method attrid=models;
  title1 "Figure 3 - Treatment Effect Model SE (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Average Model SE" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;


*** FIG 4 TREATMENT EFFECT CI COVERAGE - ADJUSTED ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="FIG4";
ods html path=out gpath=out file="fig4.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_cic_mean2 / group = method attrid=models;
  title1 "Figure 4 - Treatment Effect 95% CI Coverage (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "95% CI Coverage" values=(0.85 to 1.00 by 0.02) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;

/*  */
/* *** FIG 5 TREATMENT EFFECT BASIC SUPERIORITY POWER - ADJUSTED ***; */
/* ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="FIG5"; */
/* *ods html path=out gpath=out file="fig5.html"; */
/*  */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_pwr_mean2 / group = method attrid=models; */
/*   title1 "Figure 5 - Treatment Effect Basic Superiority Power (adjusted*)"; */
/*   footnote1 "* Variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Power" values=(0.9 to 1.00 by 0.01) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */

*** FIG 5 TREATMENT EFFECT SUPER SUPERIORITY POWER - ADJUSTED ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="FIG5";
ods html path=out gpath=out file="fig5.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_pwr2_mean2 / group = method attrid=models;
  title1 "Figure 5 - Treatment Effect Super Superiority Power (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Power" values=(0.5 to 1.00 by 0.05) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;








ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="FIG6";
proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_hw_mean2 / group = method attrid=models;
  title1 "Figure X - Treatment Effect Halfwidth (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Power" values=(0.20 to 0.35 by 0.025) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;





********************************************************************;
*** CONFIRMED SUPPLIMENTARY INFORMATION FIGURES                  ***;
********************************************************************;


*** FIG 4 TREATMENT EFFECT CI COVERAGE - ADJUSTED ***;
/* ods graphics on / reset=all width = 10in height=7in outputfmt=svg; */
/* ods html path=out gpath=out file="SI-Figures.html"; */

*** SI FIG 1 TREATMENT EFFECT BIAS - DNAR SCENARIOS ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="SI-FIG1";
ods html path=out gpath=out file="SI-fig1.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (5 6 7 8))) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_bias_mean / group = method attrid=models;
  title1 "SI Figure 1 - Treatment effect bias";
  rowaxis label = "Bias" values=(-0.10 to 0.04 by 0.02) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;


*** SI FIG 2 TREATMENT EFFECT SD(POINT ESTIMATE) ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="SI-FIG2";
ods html path=out gpath=out file="SI-fig2.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (5 6 7 8))) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_est_sd / group = method attrid=models;
  title1 "SI Figure 2 - SD of treatement effect estimates";
  rowaxis label = "SD of treatment effect estimate" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;


*** SI FIG 3 TREATMENT EFFECT MODEL SE - ADJUSTED ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="SI-FIG3";
ods html path=out gpath=out file="SI-fig3.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (5 6 7 8) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_se_mean2 / group = method attrid=models;
  title1 "SI Figure 3 - Treatment Effect Model SE (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Average Model SE" values=(0.1 to 0.18 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;


*** SI FIG 4 TREATMENT EFFECT CI COVERAGE - ADJUSTED ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="SI-FIG4";
ods html path=out gpath=out file="SI-fig4.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (5 6 7 8) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_cic_mean2 / group = method attrid=models;
  title1 "SI Figure 4 - Treatment Effect 95% CI Coverage (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "95% CI Coverage" values=(0.85 to 1.00 by 0.02) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;



********************************************************************;
*** POSSIBLE SUPPLIMENTARY INFORMATION FIGURES (ADJUSTED)        ***;
********************************************************************;

*** SI FIG 5 TREATMENT EFFECT ROOT OF AVERAGE MSE - ADJUSTED ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="SI-FIG5";
ods html path=out gpath=out file="SI-fig5.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4)  )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_rmse_mean2 / group = method attrid=models;
  title1 "SI Figure 5 - Treatment Effect Root of Average MSE (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Root of Average MSE" values=(0.15 to 0.25 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;

*** SI FIG 6 TREATMENT EFFECT ROOT OF AVERAGE MSE - ADJUSTED ***;
ods graphics on / reset=all width = 10in height=7in outputfmt=svg imagename="SI-FIG6";
ods html path=out gpath=out file="SI-fig6.html";

proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (5 6 7 8) )) dattrmap=map1;
  panelby typen / rows=2 columns=2 novarname;
  series x=rdrate y=dif_policy_rmse_mean2 / group = method attrid=models;
  title1 "SI Figure 6 - Treatment Effect Root of Average MSE (adjusted*)";
  footnote1 "* Variance correction applied to MMRM2 or MMRM3";
  rowaxis label = "Root of Average MSE" values=(0.15 to 0.25 by 0.01) grid;
  colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;  
  format rdrate rdrate. typen typen.;
run;




















/*  */
/*  */
/* *** SI FIG A.3 DIFFERENCE IN TREATMENT EFFECT SD AND MODEL ESTIMATED SE - ADJUSTED ***; */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_se_bias2 / group = method attrid=models; */
/*   title1 "SI Figure A.3 - Treatment Effect Difference in Model SE and SD of Point Estimates (adjusted*)"; */
/*   footnote1 "* Variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Difference: Model SE - SD of Point Estimates" values=(-0.005 to 0.02 by 0.005) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */
/*  */
/* *** SI FIG A.4 DIFFERENCE IN TREATMENT EFFECT SD AND MODEL ESTIMATED SE - ADJUSTED ***; */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_se_bias2 / group = method attrid=models; */
/*   title1 "SI Figure A.4 - Treatment Effect Difference in Model SE and SD of Point Estimates (adjusted*)"; */
/*   footnote1 "* Variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Difference: Model SE - SD of Point Estimates" values=(-0.005 to 0.02 by 0.005) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */
/*  */
/*  */
/*  */
/* ********************************************************************; */
/* *** POSSIBLE SUPPLIMENTARY INFORMATION FIGURES (UNADJUSTED)      ***; */
/* ********************************************************************; */
/*  */
/* *** SI FIG U.1 TREATMENT EFFECT MODEL SE - UNADJUSTED ***; */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_se_mean1 / group = method attrid=models; */
/*   title1 "SI Figure U.1 - Treatment Effect Model SE (unadjusted*)"; */
/*   footnote1 "* No variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Average Model SE" values=(0.1 to 0.18 by 0.01) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */
/*  */
/* *** SI FIG U.6 TREATMENT EFFECT MODEL SE - UNADJUSTED ***; */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (5 6 7 8) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_se_mean1 / group = method attrid=models; */
/*   title1 "SI Figure U.2 - Treatment Effect Model SE (unadjusted*)"; */
/*   footnote1 "* No variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Average Model SE" values=(0.1 to 0.18 by 0.01) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */
/*  */
/* *** SI FIG U.7 TREATMENT EFFECT ROOT OF AVERAGE MSE - UNADJUSTED ***; */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_rmse_mean1 / group = method attrid=models; */
/*   title1 "SI Figure U.3 - Treatment Effect Root of Average MSE (unadjusted*)"; */
/*   footnote1 "* No variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Root of Average MSE" values=(0.15 to 0.25 by 0.01) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */
/*  */
/* *** SI FIG U.4 TREATMENT EFFECT ROOT OF AVERAGE MSE - UNADJUSTED ***; */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (5 6 7 8) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_rmse_mean1 / group = method attrid=models; */
/*   title1 "SI Figure U.4 - Treatment Effect Root of Average MSE (unadjusted*)"; */
/*   footnote1 "* No variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Root of Average MSE" values=(0.15 to 0.25 by 0.01) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */
/*  */
/* *** SI FIG U.5 DIFFERENCE IN TREATMENT EFFECT SD AND MODEL ESTIMATED SE - UNADJUSTED ***; */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (1 2 3 4) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_se_bias1 / group = method attrid=models; */
/*   title1 "SI Figure U.5 - Treatment Effect Difference in Model SE and SD of Point Estimates (unadjusted*)"; */
/*   footnote1 "* No variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Difference: Model SE - SD of Point Estimates" values=(-0.005 to 0.02 by 0.005) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */
/*  */
/* *** SI FIG U.6 DIFFERENCE IN TREATMENT EFFECT SD AND MODEL ESTIMATED SE - UNADJUSTED ***; */
/* proc sgpanel data = all_dif_plot (where = (visitn = 5 and typen in (5 6 7 8) )) dattrmap=map1; */
/*   panelby typen / rows=2 columns=2 novarname; */
/*   series x=rdrate y=dif_policy_se_bias1 / group = method attrid=models; */
/*   title1 "SI Figure U.6 - Treatment Effect Difference in Model SE and SD of Point Estimates (unadjusted*)"; */
/*   footnote1 "* No variance correction applied to MMRM2 or MMRM3"; */
/*   rowaxis label = "Difference: Model SE - SD of Point Estimates" values=(-0.005 to 0.02 by 0.005) grid; */
/*   colaxis label = "Missingness scenario" values=(1 to 6 by 1) grid;   */
/*   format rdrate rdrate. typen typen.; */
/* run; */






















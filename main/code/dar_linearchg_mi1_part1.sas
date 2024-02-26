/*******************************************************************
| Name       : dar_linearchg_mi1_part1.sas
| Purpose    : Code to fit MI analysis 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 14MAR22 
|-------------------------------------------------------------------
| Notes:
|
| Reads in data and fits MI1 model to the data to create MI datasets
| Runs an ANCOVA on each RDRATE, SIM_RUN, IMPUTATION AND VISITN 
| fitting CHANGE as dependent variable with GROUPN and BASELINE_VAR 
| terms. LSM and DIF results then combined using Rubins rules.
| 
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%let scenario = dar_linearchg;
%let part     = 1;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";


********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

data ds_wide;
  set data.&scenario._data;
  where sim_run le 2500;
run;


********************************************************************;
*** FIT MI1 - NO ON/OFF INDICATORS                               ***;
********************************************************************;

%macro run_mi1(rdrate=, seed=, nimp=);

%ods_off(notes=N);
proc mi data    = ds_wide (where = (rdrate = &rdrate.))
        out     = mi1 (rename = _imputation_ = imputation)
        nimpute = &nimp.
        seed    = &seed.;
  by rdrate sim_run;
  var groupn baseline_var y1 y2 y3 y4 y5;
  monotone reg (y1 = groupn baseline_var);
  monotone reg (y2 = groupn y1 baseline_var);
  monotone reg (y3 = groupn y2 y1 baseline_var);
  monotone reg (y4 = groupn y3 y2 y1 baseline_var);
  monotone reg (y5 = groupn y4 y3 y2 y1 baseline_var);
run;
%ods_on();


********************************************************************;
*** TRANSPOSE MI DATA BACK INTO LONG FORMAT FOR MIXED            ***;
********************************************************************;

data mi_long;
  set mi1;
  by rdrate sim_run imputation groupn;

  array y[5] y1-y5;   *** ARRAY TO HOLD RESPONSES ***;
  array d[5] d1-d5;   *** ARRAY TO HOLD DISCONTINUATION INDICATORS ***;
  array w[5] w1-w5;   *** ARRAY TO HOLD WITHDRAWAL INDICATORS ***;

  do j = 1 to 5;
    visitn   = j;
    response = y[j];
    disc     = d[j];
    with     = w[j];
	change   = response - baseline_var;
    output;
  end;

  keep rdrate sim_run subjid groupn group imputation visitn baseline_var 
       response change disc disctime with withtime;

run;


********************************************************************;
*** RUN ANCOVA ON EACH RDRATE AND SIMULATION AND IMPUTATION      ***;
********************************************************************;

proc sort data = mi_long;
  by rdrate sim_run imputation visitn;
run;

%ods_off(notes=N);
proc mixed data = mi_long;
  by rdrate sim_run imputation visitn;
  class groupn;
  model change = groupn baseline_var / noint;
  lsmeans groupn / diff=all cl alpha=0.05;
  ods output lsmeans = lsm;
  ods output diffs   = dif;
run;
%ods_on();

proc datasets lib = work nolist;
  delete mi: ;                *** REMOVE SUBJECT LEVEL DATA ***;
quit;


********************************************************************;
*** COMBINE LSM RESULTS USING RUBINS RULES                       ***;
********************************************************************;

proc sort data = lsm;
  by rdrate sim_run groupn visitn imputation;
run;

%ods_off();
proc mianalyze data = lsm;
  by rdrate sim_run groupn visitn;
  modeleffects estimate;
  stderr stderr;
  ods output parameterestimates = lsm_mia;
run;
%ods_on();

data lsm_policy_&rdrate.;
  set lsm_mia;
  miseed = &seed.;
  lsm_policy_nimp  = nimpute;
  lsm_policy_est   = estimate;
  lsm_policy_var   = stderr**2;
  lsm_policy_se    = stderr;
  lsm_policy_lower = ifn( lclmean ne ., lclmean, estimate - probit(0.975)*stderr);  *** IF NO MISSING DATA LCLMEAN IS MISSING ***;
  lsm_policy_upper = ifn( uclmean ne ., uclmean, estimate + probit(0.975)*stderr);  *** IF NO MISSING DATA UCLMEAN IS MISSING ***;
  keep rdrate sim_run visitn groupn lsm_: miseed;
run;


********************************************************************;
*** COMBINE DIF RESULTS USING RUBINS RULES                       ***;
********************************************************************;

proc sort data = dif;
  by rdrate sim_run groupn _groupn visitn imputation;
run;

%ods_off();
proc mianalyze data = dif;
  by rdrate sim_run groupn _groupn visitn;
  modeleffects estimate;
  stderr stderr;
  ods output parameterestimates = dif_mia;
run;
%ods_on();

data dif_policy_&rdrate.;
  set dif_mia;
  miseed = &seed.;
  dif_policy_est   = estimate;
  dif_policy_var   = stderr**2;
  dif_policy_se    = stderr;
  dif_policy_lower = ifn( lclmean ne ., lclmean, estimate - probit(0.975)*stderr);  *** IF NO MISSING DATA LCLMEAN IS MISSING ***;
  dif_policy_upper = ifn( uclmean ne ., uclmean, estimate + probit(0.975)*stderr);  *** IF NO MISSING DATA UCLMEAN IS MISSING ***;
  _visitn           = visitn;
  keep rdrate sim_run groupn _groupn visitn _visitn dif_: miseed;
run;

%mend run_mi1;


********************************************************************;
*** CALL ONCE FOR EACH RDRATE                                    ***;
********************************************************************;

%run_mi1(rdrate=1, seed=1111, nimp=25);
%run_mi1(rdrate=2, seed=2222, nimp=25);
%run_mi1(rdrate=3, seed=3333, nimp=25);
%run_mi1(rdrate=4, seed=4444, nimp=25);
%run_mi1(rdrate=5, seed=5555, nimp=25);
%run_mi1(rdrate=6, seed=6666, nimp=25);


********************************************************************;
*** STACK ALL THE RESULTS BACK TOGETHER                          ***;
********************************************************************;

data lsm_policy;
  set lsm_policy_1
      lsm_policy_2
      lsm_policy_3
      lsm_policy_4
      lsm_policy_5
      lsm_policy_6;
run;

data dif_policy;
  set dif_policy_1
      dif_policy_2
      dif_policy_3
      dif_policy_4
      dif_policy_5
      dif_policy_6;
run;

********************************************************************;
*** ORGANISE RESULTS AND OUTPUT TO PERM LOCATION                 ***;
********************************************************************;

proc sort data = lsm_policy (keep = rdrate sim_run groupn visitn lsm_policy_: miseed)
          out  = results.&scenario._mi1_lsm_part&part.;
  by rdrate sim_run groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn dif_policy_: miseed)
          out  = results.&scenario._mi1_dif_part&part.;
  by rdrate sim_run visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY AND REMOVE TEMP MI DATA FROM RESULTS    ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: mi: lsm: dif: ;
quit;





















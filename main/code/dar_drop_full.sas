/*******************************************************************
| Name       : dar_drop_full.sas
| Purpose    : Code to fit MMRM to full sim study data 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 15SEP21 
|-------------------------------------------------------------------
| Notes:
|
| S1: Reads in a the simulated data.
|
| S2. Analysis: This is an MMRM fitted to CHANGE_FULL 
|     but with no covariate to make distinctions between on and off 
|     treatment values. As there is no missing data this naive model
|     is not biased. 
| 
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%let scenario = dar_drop;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";


********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

data ds;
  set data.&scenario.;
  where sim_run le 10000;
run;


********************************************************************;
*** ANALYSIS: MMRM ON FULL DATA                                  ***;
********************************************************************;

%ods_off();
proc mixed data = ds noclprint;
  by rdrate sim_run;
  class groupn visitn subjid;
  model change_full = groupn*visitn baseline_var*visitn / noint;
  repeated visitn / subject=subjid type=un;
  lsmeans groupn*visitn / diff=all cl alpha=0.05;
  ods output lsmeans = lsm;
  ods output diffs = dif (where = (visitn = _visitn));
run;
%ods_on();


********************************************************************;
*** SECTION 3: ORGANISE RESULTS AND OUTPUT TO PERM LOCATION      ***;
********************************************************************;

data results.&scenario._full_lsm;
  set lsm;
  lsm_policy_est   = estimate;
  lsm_policy_se    = stderr;
  lsm_policy_lower = lower;
  lsm_policy_upper = upper;
  keep rdrate sim_run groupn visitn lsm_policy_:;
run;

data results.&scenario._full_dif;
  set dif;
  dif_policy_est   = estimate;
  dif_policy_se    = stderr;
  dif_policy_lower = lower;
  dif_policy_upper = upper;
  keep rdrate sim_run groupn _groupn visitn _visitn dif_policy_:;
run;


********************************************************************;
*** SECTION 3: TIDY UP WORK LIBRARY                              ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: lsm: dif: ;
quit;

/*******************************************************************
| Name       : dnar_linearchg_true.sas
| Purpose    : Code to fit MMRM to all data to act as the truth 
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
|     is not biased. This is fitted to all the simulations for each 
|     RDRATE (200000 subjects per arm) to act as the truth.
| 
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%let scenario = dnar_linearchg;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";

********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

data ds;
  set data.&scenario._true;
  where sim_run le 5000;
  all = 1;
run;


********************************************************************;
*** MMRM FOR EACH RDRATE                                         ***;
********************************************************************;

%ods_off();
proc mixed data = ds noclprint;
  by all;
  class groupn visitn usubjid;
  model change_full = groupn*visitn baseline_var*visitn / noint;
  repeated visitn / subject=usubjid type=un;
  lsmeans groupn*visitn / diff=all cl alpha=0.05;
  ods output lsmeans = lsm;
  ods output diffs = dif (where = (visitn = _visitn));
run;
%ods_on();


********************************************************************;
*** ORGANISE RESULTS AND OUTPUT TO PERM LOCATION                 ***;
********************************************************************;

data results.&scenario._true_lsm;
  set lsm;
  lsm_true_est   = estimate;
  lsm_true_se    = stderr;
  lsm_true_lower = lower;
  lsm_true_upper = upper;
  keep all groupn visitn df lsm_true_:;
run;

data results.&scenario._true_dif;
  set dif;
  dif_true_est   = estimate;
  dif_true_se    = stderr;
  dif_true_lower = lower;
  dif_true_upper = upper;
  keep all groupn visitn df dif_true_:;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY                                         ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: lsm: dif: ;
quit;

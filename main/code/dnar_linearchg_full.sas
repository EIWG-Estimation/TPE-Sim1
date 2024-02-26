/*******************************************************************
| Name       : dnar_linearchg_full.sas
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

%let scenario = dnar_linearchg;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";


********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

data ds_wide;
  set data.&scenario._data;
  where sim_run le 5000;
run;

********************************************************************;
*** CREATE LONG FORMAT DATA                                      ***;
********************************************************************;

data ds_long;
  set ds_wide;
  by rdrate sim_run;
  
  array yf[5] yf1-yf5;  *** ARRAY TO HOLD FULL RESPONSES ***;

  array y[5] y1-y5;     *** ARRAY TO HOLD RESPONSES ***;
  array d[5] d1-d5;     *** ARRAY TO HOLD DISCONTINUATION INDICATORS ***;
  array p[5] p1-p5;     *** ARRAY TO HOLD PATTERNS AT EACH VISIT ***;
  array w[5] w1-w5;     *** ARRAY TO HOLD WITHDRAWAL INDICATORS ***;
  
  array disc_code[5] disc1-disc5;     
  array pat_code[5] pat1-pat5;

  do j = 1 to 5;
    visitn        = j;
    response_full = yf[j];
    response      = y[j];
    disc          = d[j];
    pattern       = p[j];
    with          = w[j];
    disccode      = disc_code[j];
    patcode       = pat_code[j];
    if response_full ne . then change_full = response_full - baseline_var;
    else change_full = .;
    if response ne . then change = response - baseline_var;
    else change = .;
    output;
  end;

  keep rdrate sim_run subjid groupn group visitn patmaxn patmax
       baseline_var response_full response change_full change 
       disc disccode disctime pattern patcode with withtime 
       p_i1-p_i5 d_i1-d_i5 e_i1-e_i5;

run;


********************************************************************;
*** ANALYSIS: MMRM ON FULL DATA                                  ***;
********************************************************************;

%ods_off();
proc mixed data = ds_long noclprint;
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
  lsm_policy_var   = stderr**2;
  lsm_policy_lower = lower;
  lsm_policy_upper = upper;
  keep rdrate sim_run groupn visitn lsm_policy_:;
run;

data results.&scenario._full_dif;
  set dif;
  dif_policy_est   = estimate;
  dif_policy_se    = stderr;
  dif_policy_var   = stderr**2;
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

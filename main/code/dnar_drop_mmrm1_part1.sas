/*******************************************************************
| Name       : dnar_drop_mmrm1_part1.sas
| Purpose    : Code to fit naive MMRM to sim study data 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 15SEP21 
|-------------------------------------------------------------------
| Notes:
| 
| This is a basic "naive" MMRM fitted to CHANGE but with no 
| covariate to make distinctions between on and off treatment 
| values. This is biased.
| 
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                          ***;
********************************************************************;

%let scenario = dnar_drop;
%let part     = 1;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";

********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

data ds_wide;
  set data.&scenario._data;
  where 0 lt sim_run le 2500;
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
*** FIT MMRM1 TO DATA WITH CORRECT BASELINE MEAN                 ***;
********************************************************************;

%macro mmrm1(rdrate=);

  %ods_off();

  *** COUNT SIMS ***;
  data ds1;
    set ds_long (where = (rdrate = &rdrate.)) end = eof;
    retain sim_count 0;
    by rdrate sim_run;
    if first.sim_run then sim_count + 1;
	if eof then call symput("nsims", put(sim_count,best.));
  run;

  %do ii = 1 %to &nsims.;

  *** SQL TO CREATE MACRO VARIABLE FOR BASELINE MEAN ***;
  proc sql noprint;
    select mean(baseline_var) into :mbl
    from ds1
    where sim_count = &ii. and visitn = 1;
  quit;

  *** MMRM1 WITH LSM AT CORRECT BASELINE VALUE ***;
  proc mixed data = ds1 (where = (sim_count = &ii.));
    by rdrate sim_run;
    class visitn groupn subjid;
    model change = visitn*groupn visitn*baseline_var / noint ddfm=kr;
    repeated visitn / subject=subjid type=un;
    lsmeans visitn*groupn / diff=all at baseline_var = &mbl. cl;
    ods output lsmeans = lsm_sc&ii.;
    ods output diffs   = dif_sc&ii. (where = (visitn = _visitn));
  run;

  *** STACK LSM RESULTS FOR EACH SIM_RUN ***; 
  data lsm_r&rdrate.;
    set %if &ii. ne 1 %then %do; lsm_r&rdrate. %end;
        lsm_sc&ii.;
    by rdrate sim_run visitn groupn;
  run;

  *** STACK DIF RESULTS FOR EACH SIM_RUN ***; 
  data dif_r&rdrate.;
    set %if &ii. ne 1 %then %do; dif_r&rdrate. %end;
        dif_sc&ii.;
    by rdrate sim_run visitn groupn;
  run;

  %end;
  %ods_on();

%mend;

%mmrm1(rdrate=1);
%mmrm1(rdrate=2);
%mmrm1(rdrate=3);
%mmrm1(rdrate=4);
%mmrm1(rdrate=5);
%mmrm1(rdrate=6);


********************************************************************;
*** STACK RESULTS AND CLEAN UP WORK AREA                         ***;
********************************************************************;

data lsm_policy;
  set lsm_r1
      lsm_r2
      lsm_r3
      lsm_r4
      lsm_r5
      lsm_r6;
  by rdrate sim_run visitn groupn;

  lsm_policy_est   = estimate;
  lsm_policy_var   = stderr**2;
  lsm_policy_se    = stderr;
  lsm_policy_lower = lower;
  lsm_policy_upper = upper;

  keep rdrate sim_run visitn groupn lsm_:;

run; 

data dif_policy;
  set dif_r1
      dif_r2
      dif_r3
      dif_r4
      dif_r5
      dif_r6;
  by rdrate sim_run visitn groupn;

  dif_policy_est   = estimate;
  dif_policy_var   = stderr**2;
  dif_policy_se    = stderr;
  dif_policy_lower = lower;
  dif_policy_upper = upper;

  keep rdrate sim_run visitn _visitn groupn _groupn dif_:;

run; 


********************************************************************;
*** ORGANISE RESULTS AND OUTPUT TO PERM LOCATION                 ***;
********************************************************************;

proc sort data = lsm_policy (keep = rdrate sim_run groupn visitn lsm_policy_:)
          out  = results.&scenario._mmrm1_lsm_part&part.;
  by rdrate sim_run groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn dif_policy_:)
          out  = results.&scenario._mmrm1_dif_part&part.;
  by rdrate sim_run visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY                                         ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: lsm: dif: ;
quit;


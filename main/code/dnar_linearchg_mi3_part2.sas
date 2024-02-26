/*******************************************************************
| Name       : dnar_linearchg_mi3_part2.sas
| Purpose    : Code to fit MI analysis 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 15APR22 
|-------------------------------------------------------------------
| Notes: MI3
|
| 
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%let scenario = dnar_linearchg;
%let part     = 2;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";


********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

data ds_wide;
  set data.&scenario._data;
  where 2500 lt sim_run le 5000;
run;


********************************************************************;
*** FIT MI3 MODEL                                                ***;
********************************************************************;

%macro run_mi3(rdrate=, nimp=, seed=);

%let seed1 = %eval(&seed. + 100);
%let seed2 = %eval(&seed. + 200);
%let seed3 = %eval(&seed. + 300);
%let seed4 = %eval(&seed. + 400);
%let seed5 = %eval(&seed. + 500);


********************************************************************;
*** FIT MI3 USING STANDARDIZED RESIDUALS                         ***;
********************************************************************;

%ods_off(notes=Y);

*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi3;
  set ds_wide (where = (rdrate = &rdrate.));
  dummy = y5;  
run;

*** IMPUTE Y1 DIRECTLY USING BASELINE AS ALL ON TREATMENT ***;
proc mi data    = mi3
        out     = mi3_y1 (rename=(_imputation_ = imputation) drop = dummy)
        nimpute = &nimp.
        seed    = &seed1.;
  by rdrate sim_run patmaxn patmax;
  class groupn pat1;
  var groupn pat1 baseline_var y1 dummy;
  monotone reg (y1 = groupn pat1 groupn*pat1 baseline_var);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi3_y1;
  by rdrate sim_run imputation;
  class groupn pat1 patmaxn patmax;
  model y1 = groupn pat1 groupn*pat1 baseline_var / noint residuals
  outpm = mi3_r1 (rename=(resid=r1) drop = studentresid pearsonresid pred stderrpred df alpha lower upper) ;
run;

*** TIDY UP INTERIM DATA ***;
proc datasets lib = work nolist;
  delete mi3_y1:;
quit;


*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi3;
  set mi3_r1;
  dummy = y5;  
run;

*** IMPUTE Y2 USING STANDARDIZED RESIDUAL R1 ***;
proc mi data    = mi3 
        out     = mi3_y2 (drop = dummy)
        nimpute = 1
        seed    = &seed2.;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn pat2;
  var groupn pat2 baseline_var r1 y2 dummy;
  monotone reg (y2 = groupn pat2 groupn*pat2 baseline_var r1);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi3_y2;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn pat2;
  model y2 = groupn pat2 groupn*pat2 baseline_var r1 / noint residuals
  outpm = mi3_r2 (rename=(resid=r2) drop = studentresid pearsonresid pred stderrpred df alpha lower upper) ;
run;

*** TIDY UP INTERIM DATA ***;
proc datasets lib = work nolist;
  delete mi3_r1 mi3_y2;
quit;


*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi3;
  set mi3_r2;
  dummy = y5;  
run;

*** IMPUTE Y3 USING STANDARDIZED RESIDUALS R1 R2 ***;
proc mi data    = mi3
        out     = mi3_y3 (drop = dummy)
        nimpute = 1
        seed    = &seed3.;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn pat3;
  var groupn pat3 baseline_var r1 r2 y3 dummy;
  monotone reg (y3 = groupn pat3 groupn*pat3 baseline_var r1 r2);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi3_y3;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn pat3;
  model y3 = groupn pat3 groupn*pat3 baseline_var r1 r2 / noint residuals
  outpm = mi3_r3 (rename=(resid=r3) drop = studentresid pearsonresid pred stderrpred df alpha lower upper) ;
run;

*** TIDY UP INTERIM DATA ***;
proc datasets lib = work nolist;
  delete mi3_r2 mi3_y3;
quit;


*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi3;
  set mi3_r3;
  dummy = y5;  
run;

*** IMPUTE Y4 USING STANDARDIZED RESIDUALS R1 R2 R3 ***;
proc mi data    = mi3
        out     = mi3_y4 (drop = dummy)
        nimpute = 1
        seed    = &seed4.;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn pat4;
  var groupn pat4 baseline_var r1 r2 r3 y4 dummy;
  monotone reg (y4 = groupn pat4 groupn*pat4 baseline_var r1 r2 r3);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi3_y4;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn pat4;
  model y4 = groupn pat4 groupn*pat4 baseline_var r1 r2 r3 / noint residuals
  outpm = mi3_r4 (rename=(resid=r4) drop = studentresid pearsonresid pred stderrpred df alpha lower upper) ;
run;

*** TIDY UP INTERIM DATA ***;
proc datasets lib = work nolist;
  delete mi3_r3 mi3_y4;
quit;


*** IMPUTE Y5 USING STANDARDIZED RESIDUALS R1 R2 R3 R4 ***;
proc mi data    = mi3_r4 
        out     = mi3
        nimpute = 1
        seed    = &seed5.;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn pat5;
  var groupn pat5 baseline_var r1 r2 r3 r4 y5;
  monotone reg (y5 = groupn pat5 groupn*pat5 baseline_var r1 r2 r3 r4);
run;

proc datasets lib = work nolist;
  delete mi3_r4;
quit;

%ods_on();


********************************************************************;
*** TRANSPOSE MI DATA BACK INTO LONG FORMAT FOR MIXED            ***;
********************************************************************;

data mi3_long;
  set mi3;
  by rdrate sim_run patmaxn patmax imputation;

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
       response change disc disctime with withtime patmaxn patmax;

run;

proc datasets lib = work nolist;
  delete mi3;
quit;


********************************************************************;
*** RUN ANCOVA ON EACH SIMULATION AND IMPUTATION                 ***;
********************************************************************;

proc sort data = mi3_long;
  by rdrate sim_run patmaxn patmax visitn imputation;
run;

%ods_off(notes=N);
proc mixed data = mi3_long;
  by rdrate sim_run patmaxn patmax visitn imputation;
  class groupn;
  model change = groupn baseline_var / noint;
  lsmeans groupn / diff=all cl alpha=0.05;
  ods output lsmeans = lsm;
  ods output diffs   = dif;
run;
%ods_on();


********************************************************************;
*** COMBINE LSM RESULTS USING RUBINS RULES                       ***;
********************************************************************;

proc sort data = lsm;
  by rdrate sim_run patmaxn patmax groupn visitn imputation;
run;

%ods_off();
proc mianalyze data = lsm;
  by rdrate sim_run patmaxn patmax groupn visitn;
  modeleffects estimate;
  stderr stderr;
  ods output parameterestimates = lsm_mia;
run;
%ods_on();

data lsm_policy_&rdrate.;
  set lsm_mia;

  lsm_policy_nimp  = nimpute;
  lsm_policy_est   = estimate;
  lsm_policy_var   = stderr**2;
  lsm_policy_se    = stderr;
  lsm_policy_lower = ifn( lclmean ne ., lclmean, estimate - probit(0.975)*stderr);  *** IF NO MISSING DATA LCLMEAN IS MISSING ***;
  lsm_policy_upper = ifn( uclmean ne ., uclmean, estimate + probit(0.975)*stderr);  *** IF NO MISSING DATA UCLMEAN IS MISSING ***;

  keep rdrate sim_run patmaxn patmax visitn groupn lsm_:;

run;


********************************************************************;
*** COMBINE DIF RESULTS USING RUBINS RULES                       ***;
********************************************************************;

proc sort data = dif;
  by rdrate sim_run patmaxn patmax groupn _groupn visitn imputation;
run;

%ods_off();
proc mianalyze data = dif;
  by rdrate sim_run patmaxn patmax groupn _groupn visitn;
  modeleffects estimate;
  stderr stderr;
  ods output parameterestimates = dif_mia;
run;
%ods_on();

data dif_policy_&rdrate.;
  set dif_mia;

  dif_policy_est   = estimate;
  dif_policy_var   = stderr**2;
  dif_policy_se    = stderr;
  dif_policy_lower = ifn( lclmean ne ., lclmean, estimate - probit(0.975)*stderr);  *** IF NO MISSING DATA LCLMEAN IS MISSING ***;
  dif_policy_upper = ifn( uclmean ne ., uclmean, estimate + probit(0.975)*stderr);  *** IF NO MISSING DATA UCLMEAN IS MISSING ***;
  _visitn           = visitn;

  keep rdrate sim_run patmaxn patmax groupn _groupn visitn _visitn dif_:;

run;

%mend run_mi3;


********************************************************************;
*** CALL MI3 AND STACK RESULTS                                   ***;
********************************************************************;

%run_mi3(rdrate=1, seed=1111, nimp=25);
%run_mi3(rdrate=2, seed=2222, nimp=25);
%run_mi3(rdrate=3, seed=3333, nimp=25);
%run_mi3(rdrate=4, seed=4444, nimp=25);
%run_mi3(rdrate=5, seed=5555, nimp=25);
%run_mi3(rdrate=6, seed=6666, nimp=25);   



data lsm_policy;
  set lsm_policy_1
      lsm_policy_2
      lsm_policy_3
      lsm_policy_4
      lsm_policy_5
      lsm_policy_6  
      ;
  by rdrate sim_run patmaxn patmax groupn visitn;
  
  fitn = 6 - patmaxn + 1;
  length fit $50;
  select(fitn);
    when (1) fit = "6-pattern MI (MI3)";
    when (2) fit = "5-pattern MI";
    when (3) fit = "4-pattern MI";
    when (4) fit = "3-pattern MI";
    when (5) fit = "2-pattern MI (MI2)";
    when (6) fit = "1-pattern MI (MI1)";
    otherwise;
  end;
   
run;

data dif_policy;
  set dif_policy_1
      dif_policy_2
      dif_policy_3
      dif_policy_4
      dif_policy_5
      dif_policy_6
      ;
  by rdrate sim_run patmaxn patmax visitn;

  fitn = 6 - patmaxn + 1;
  length fit $50;
  select(fitn);
    when (1) fit = "6-pattern MI (MI3)";
    when (2) fit = "5-pattern MI";
    when (3) fit = "4-pattern MI";
    when (4) fit = "3-pattern MI";
    when (5) fit = "2-pattern MI (MI2)";
    when (6) fit = "1-pattern MI (MI1)";
    otherwise;
  end;

run;

proc datasets lib = work nolist;
  delete lsm_policy_: dif_policy_:;
quit;


********************************************************************;
*** ORGANISE RESULTS AND OUTPUT TO PERM LOCATION                 ***;
********************************************************************;

proc sort data = lsm_policy (keep = rdrate sim_run groupn visitn patmaxn patmax fitn fit lsm_policy_: )
          out  = results.&scenario._mi3_lsm_part&part.;
  by rdrate sim_run groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn patmaxn patmax fitn fit dif_policy_: )
          out  = results.&scenario._mi3_dif_part&part.;
  by rdrate sim_run visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY AND REMOVE TEMP MI DATA FROM RESULTS    ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: mi: lsm: dif: ;
quit;














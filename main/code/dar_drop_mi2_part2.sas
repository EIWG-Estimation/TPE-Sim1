/*******************************************************************
| Name       : dar_drop_mi2_part2.sas
| Purpose    : Code to fit MI analysis 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 14MAR22 
|-------------------------------------------------------------------
| Notes:
|
| Reads in data and works out if it can fit MI2 or only MI1 model 
| to the data.
| Fits either MI2 or MI1 model and uses it to create MI datasets.
| Runs an ANCOVA on each RDRATE, SIM_RUN, IMPUTATION AND VISITN 
| fitting CHANGE as dependent variable with GROUPN and BASELINE_VAR 
| terms. LSM and DIF results then combined using Rubins rules.
| 
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%let scenario = dar_drop;
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
*** FIT MI2 MODEL                                                ***;
********************************************************************;

%macro run_mi2(rdrate=, nimp=, seed=);

%let seed1 = %eval(&seed. + 100);
%let seed2 = %eval(&seed. + 200);
%let seed3 = %eval(&seed. + 300);
%let seed4 = %eval(&seed. + 400);
%let seed5 = %eval(&seed. + 500);


********************************************************************;
*** FIT MI2 USING STANDARDIZED RESIDUALS                         ***;
********************************************************************;

%ods_off(notes=Y);


*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi2;
  set ds_wide (where = (rdrate = &rdrate.));
  dummy = y5;  
run;

*** IMPUTE Y1 DIRECTLY USING BASELINE AS ALL ON TREATMENT ***;
proc mi data    = mi2
        out     = mi2_y1 (rename=(_imputation_ = imputation) drop = dummy)
        nimpute = &nimp.
        seed    = &seed1.;
  by rdrate sim_run patmaxn patmax;
  class groupn disc1;
  var groupn disc1 baseline_var y1 dummy;
  monotone reg (y1 = groupn disc1 groupn*disc1 baseline_var);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi2_y1;
  by rdrate sim_run imputation;
  class groupn disc1 patmaxn patmax;
  model y1 = groupn disc1 groupn*disc1 baseline_var / noint residuals
  outpm = mi2_r1 (rename=(resid=r1) drop = studentresid pearsonresid pred stderrpred df alpha lower upper) ;
run;

*** TIDY UP INTERIM DATA ***;
proc datasets lib = work nolist;
  delete mi2_y1:;
quit;


*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi2;
  set mi2_r1;
  dummy = y5;  
run;

*** IMPUTE Y2 USING STANDARDIZED RESIDUAL R1 ***;
proc mi data    = mi2 
        out     = mi2_y2 (drop = dummy)
        nimpute = 1
        seed    = &seed2.;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn disc2;
  var groupn disc2 baseline_var r1 y2 dummy;
  monotone reg (y2 = groupn disc2 groupn*disc2 baseline_var r1);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi2_y2;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn disc2;
  model y2 = groupn disc2 groupn*disc2 baseline_var r1 / noint residuals
  outpm = mi2_r2 (rename=(resid=r2) drop = studentresid pearsonresid pred stderrpred df alpha lower upper) ;
run;

*** TIDY UP INTERIM DATA ***;
proc datasets lib = work nolist;
  delete mi2_r1 mi2_y2;
quit;


*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi2;
  set mi2_r2;
  dummy = y5;  
run;

*** IMPUTE Y3 USING STANDARDIZED RESIDUALS R1 R2 ***;
proc mi data    = mi2
        out     = mi2_y3 (drop = dummy)
        nimpute = 1
        seed    = &seed3.;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn disc3;
  var groupn disc3 baseline_var r1 r2 y3 dummy;
  monotone reg (y3 = groupn disc3 groupn*disc3 baseline_var r1 r2);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi2_y3;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn disc3;
  model y3 = groupn disc3 groupn*disc3 baseline_var r1 r2 / noint residuals
  outpm = mi2_r3 (rename=(resid=r3) drop = studentresid pearsonresid pred stderrpred df alpha lower upper) ;
run;

*** TIDY UP INTERIM DATA ***;
proc datasets lib = work nolist;
  delete mi2_r2 mi2_y3;
quit;


*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi2;
  set mi2_r3;
  dummy = y5;  
run;

*** IMPUTE Y4 USING STANDARDIZED RESIDUALS R1 R2 R3 ***;
proc mi data    = mi2
        out     = mi2_y4 (drop = dummy)
        nimpute = 1
        seed    = &seed4.;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn disc4;
  var groupn disc4 baseline_var r1 r2 r3 y4 dummy;
  monotone reg (y4 = groupn disc4 groupn*disc4 baseline_var r1 r2 r3);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi2_y4;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn disc4;
  model y4 = groupn disc4 groupn*disc4 baseline_var r1 r2 r3 / noint residuals
  outpm = mi2_r4 (rename=(resid=r4) drop = studentresid pearsonresid pred stderrpred df alpha lower upper) ;
run;

*** TIDY UP INTERIM DATA ***;
proc datasets lib = work nolist;
  delete mi2_r3 mi2_y4;
quit;


*** IMPUTE Y5 USING STANDARDIZED RESIDUALS R1 R2 R3 R4 ***;
proc mi data    = mi2_r4 
        out     = mi2
        nimpute = 1
        seed    = &seed5.;
  by rdrate sim_run patmaxn patmax imputation;
  class groupn disc5;
  var groupn disc5 baseline_var r1 r2 r3 r4 y5;
  monotone reg (y5 = groupn disc5 groupn*disc5 baseline_var r1 r2 r3 r4);
run;

proc datasets lib = work nolist;
  delete mi2_r4;
quit;

%ods_on();


********************************************************************;
*** TRANSPOSE MI DATA BACK INTO LONG FORMAT FOR MIXED            ***;
********************************************************************;

data mi2_long;
  set mi2;
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
  delete mi2;
quit;


********************************************************************;
*** RUN ANCOVA ON EACH SIMULATION AND IMPUTATION                 ***;
********************************************************************;

proc sort data = mi2_long;
  by rdrate sim_run patmaxn patmax visitn imputation;
run;

%ods_off(notes=N);
proc mixed data = mi2_long;
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

%mend run_mi2;


********************************************************************;
*** CALL MI2 AND STACK RESULTS                                   ***;
********************************************************************;

%run_mi2(rdrate=1, seed=1111, nimp=25);
%run_mi2(rdrate=2, seed=2222, nimp=25);
%run_mi2(rdrate=3, seed=3333, nimp=25);
%run_mi2(rdrate=4, seed=4444, nimp=25);
%run_mi2(rdrate=5, seed=5555, nimp=25);
%run_mi2(rdrate=6, seed=6666, nimp=25);   


data lsm_policy;
  set lsm_policy_1
      lsm_policy_2
      lsm_policy_3
      lsm_policy_4
      lsm_policy_5
      lsm_policy_6  
      ;
  by rdrate sim_run patmaxn patmax groupn visitn;
  
  fitn = ifn(patmaxn > 1, 1, 2);
  length fit $50;
  select(fitn);
    when (1) fit = "2-pattern MI (MI2)";
    when (2) fit = "1-pattern MI (MI1)";
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
 
  fitn = ifn(patmaxn > 1, 1, 2);
  length fit $50;
  select(fitn);
    when (1) fit = "2-pattern MI (MI2)";
    when (2) fit = "1-pattern MI (MI1)";
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
          out  = results.&scenario._mi2_lsm_part&part.;
  by rdrate sim_run groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn patmaxn patmax fitn fit dif_policy_: )
          out  = results.&scenario._mi2_dif_part&part.;
  by rdrate sim_run visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY AND REMOVE TEMP MI DATA FROM RESULTS    ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: mi: lsm: dif: ;
quit;


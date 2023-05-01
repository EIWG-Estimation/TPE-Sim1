/*******************************************************************
| Name       : dnar_linearchg_mi2_part1.sas
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

%let scenario = dnar_linearchg;
%let part     = 1;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";

%macro run_mi(rdrate=, nimp=, seed=);

%let seed1 = %eval(&seed. + 1);
%let seed2 = %eval(&seed. + 2);
%let seed3 = %eval(&seed. + 3);
%let seed4 = %eval(&seed. + 4);
%let seed5 = %eval(&seed. + 5);

********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

proc sql noprint;
  create table ds as
    select *, 1-ontrt as disc 
    from data.&scenario.
    where rdrate = &rdrate. and (0 lt sim_run le 5000) 
    order by rdrate, sim_run, groupn, subjid, visitn, disc;
quit;


********************************************************************;
*** TRANSPOSE TO WIDE FORMAT                                     ***;
********************************************************************;

data ds1;
  set ds;
  by rdrate sim_run groupn subjid;
  retain row 0 y1-y5 d1-d5 w1-w5 .;   *** RETAIN VARIABLES OVER ROWS ***;

  array y[5] y1-y5;   *** ARRAY TO HOLD RESPONSES ***;
  array d[5] d1-d5;   *** ARRAY TO HOLD DISCONTINUATION INDICATORS ***;
  array w[5] w1-w5;   *** ARRAY TO HOLD WITHDRAWAL INDICATORS ***;

  if first.subjid then do;   *** FOR EACH SUBJECT RESET THE ARRAYS ***;
    row = 0;
	do j = 1 to 5;
     y[j] = .;
	 d[j] = .;
	 w[j] = .;
	end;
  end;

  row    = row + 1;     *** COUNT THE ROW ***;
  y[row] = response;    
  d[row] = disc;        *** DISCONTINUED FLAG ***;
  w[row] = is_missing;  *** WITHDRAWAL FLAG ***;

  if last.subjid then do;   *** ONLY OUTPUT FINAL ROW PER PATIENT ***;
    disctime = 5 - sum(of d1-d5);
    withtime = 5 - sum(of w1-w5);
    output;
  end;
  keep rdrate sim_run subjid groupn group baseline_var disctime withtime y1-y5 d1-d5 w1-w5;

run;


********************************************************************;
*** CHECK EACH SIM FOR COMPLETE DATA AND IF IMPUTATION POSSIBLE  ***;
********************************************************************;

data ds2;
  set ds1;
  by rdrate sim_run groupn;
  retain row i1-i5 0 c1-c5 mi2_any mi2_v5 1;  

  array d[5] d1-d5;
  array w[5] w1-w5;
  array c[5] c1-c5;       *** RETAINED COMPLETE DATA FLAGS ***;
  array i[5] i1-i5;       *** RETAINED IMPUTATION FLAGS ***;

  if first.sim_run then do;
     do j = 1 to 5;
       c[j] = 1;          *** RESET THE COMPLETE DATA FLAGS TO YES ***;
	 end;
     mi2_any = 1;         *** RESET MI2 ANY MODEL FLAG AS YES ***;
     mi2_v5  = 1;         *** RESET MI2 V5 MODEL FLAG AS YES ***;
  end;

  if first.groupn then do j = 1 to 5;  
    i[j] = 0;             *** RESET BY GROUP IMPUTATION CHECK FLAGS AS NO ***;
  end;

  do j = 1 to 5;
    if w[j]=1 then c[j] = 0;            *** IF ONE SUBJECT MISS THEN DATA NOT COMPLETE ***;
    if d[j]=1 and w[j]=0 then i[j] = 1; *** IF ONE SUBJECT DISC AND DID NOT WITHDRAW IMP IS POSSIBLE ***;
  end;
 
  if last.groupn then do;                    *** CHECK EACH TREATMENT CANNOT RUN IF EITHER FAILED ***;
    if sum(of i1-i5) lt 5 then mi2_any = 0;  *** IF ANY FAILS SET MI2 ANY FLAG TO FAIL ***;
    if i5 ne 1 then mi2_v5 = 0;              *** IF V5 FAILS SET MI2 V5 FLAG TO FAIL  ***;
  end;

  if last.sim_run then output;    *** OUTPUT THE COMPLETE DATA AND MI2 FLAGS FOR EACH SIMULATION ***;
  keep rdrate sim_run c1-c5 mi2:;

run;


********************************************************************;
*** MERGE ON MODEL FITTING FLAG                                  ***;
********************************************************************;

data ds3;
  merge ds1
        ds2;
  by rdrate sim_run;
run;

proc datasets lib = work nolist;
  delete ds ds1 ds2;
quit;


********************************************************************;
*** FIT MI2 USING STANDARDIZED RESIDUALS                         ***;
********************************************************************;

%ods_off(notes=Y);

*** DUMMY VAR TO ENSURE SOME IMPUTATION ***;
data mi2;
  set ds3 (where = (mi2_any = 1));
  dummy = y5;  
run;

*** IMPUTE Y1 DIRECTLY USING BASELINE AS ALL ON TREATMENT ***;
proc mi data    = mi2
        out     = mi2_y1 (rename=(_imputation_ = imputation) drop = dummy)
        nimpute = &nimp.
        seed    = &seed1.;
  by rdrate sim_run;
  class groupn d1;
  var groupn d1 baseline_var y1 dummy;
  monotone reg (y1 = groupn d1 groupn*d1 baseline_var);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi2_y1;
  by rdrate sim_run imputation;
  class groupn d1;
  model y1 = groupn d1 groupn*d1 baseline_var / noint residuals
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
  by rdrate sim_run imputation;
  class groupn d2;
  var groupn d2 baseline_var r1 y2 dummy;
  monotone reg (y2 = groupn d2 groupn*d2 baseline_var r1);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi2_y2;
  by rdrate sim_run imputation;
  class groupn d2;
  model y2 = groupn d2 groupn*d2 baseline_var r1 / noint residuals
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
  by rdrate sim_run imputation;
  class groupn d3;
  var groupn d3 baseline_var r1 r2 y3 dummy;
  monotone reg (y3 = groupn d3 groupn*d3 baseline_var r1 r2);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi2_y3;
  by rdrate sim_run imputation;
  class groupn d3;
  model y3 = groupn d3 groupn*d3 baseline_var r1 r2 / noint residuals
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
  by rdrate sim_run imputation;
  class groupn d4;
  var groupn d4 baseline_var r1 r2 r3 y4 dummy;
  monotone reg (y4 = groupn d4 groupn*d4 baseline_var r1 r2 r3);
  monotone reg (dummy = baseline_var); *** TRICK TO KEEP MI HAPPY ***;
run;

*** GET RESIDUALS BY FITTING SAME MODEL TO IMPUTED DATA IN MIXED ***;
proc mixed data = mi2_y4;
  by rdrate sim_run imputation;
  class groupn d4;
  model y4 = groupn d4 groupn*d4 baseline_var r1 r2 r3 / noint residuals
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
  by rdrate sim_run imputation;
  class groupn d5;
  var groupn d5 baseline_var r1 r2 r3 r4 y5;
  monotone reg (y5 = groupn d5 groupn*d5 baseline_var r1 r2 r3 r4);
run;

proc datasets lib = work nolist;
  delete mi2_r4;
quit;

%ods_on();


********************************************************************;
*** FIT MI1 - NO ON/OFF INDICATORS SO NO NEED FOR RESIDUALS      ***;
********************************************************************;

%ods_off(notes=N);
proc mi data    = ds3 (where = (mi2_any = 0))
        out     = mi1 (rename = (_imputation_ = imputation))
        nimpute = &nimp.
        seed    = &seed.;
  by rdrate sim_run;
  class groupn;
  var groupn baseline_var y1 y2 y3 y4 y5;
  monotone reg (y1 = groupn baseline_var);
  monotone reg (y2 = groupn y1 baseline_var);
  monotone reg (y3 = groupn y2 y1 baseline_var);
  monotone reg (y4 = groupn y3 y2 y1 baseline_var);
  monotone reg (y5 = groupn y4 y3 y2 y1 baseline_var);
run;
%ods_on();


********************************************************************;
*** STACK ALL IMPUTATIONS AND CREATE FIT VARIABLES               ***;
********************************************************************;

data mi;
  set mi1 (in = in1 )
      mi2 (in = in2 );
  by rdrate sim_run imputation;

  length fit $50;
  if in2 then do;
    fitn = 1;
	fit  = "Full MI2";
  end;
  else if in1 then do;
    fitn = 2;
	fit  = "Reduced to MI1";
  end;

run;

proc datasets lib = work nolist;
  delete mi1 mi2;
quit;

********************************************************************;
*** TRANSPOSE MI DATA BACK INTO LONG FORMAT FOR MIXED            ***;
********************************************************************;

data mi_long;
  set mi;
  by rdrate sim_run imputation;

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
       response change disc disctime with withtime fitn fit;

run;

proc datasets lib = work nolist;
  delete mi;
quit;


********************************************************************;
*** RUN ANCOVA ON EACH RDRATE AND SIMULATION AND IMPUTATION      ***;
********************************************************************;

proc sort data = mi_long;
  by rdrate sim_run fitn fit visitn imputation;
run;

%ods_off(notes=N);
proc mixed data = mi_long;
  by rdrate sim_run fitn fit visitn imputation;
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
  by rdrate sim_run fitn fit groupn visitn imputation;
run;

%ods_off();
proc mianalyze data = lsm;
  by rdrate sim_run fitn fit groupn visitn;
  modeleffects estimate;
  stderr stderr;
  ods output parameterestimates = lsm_mia;
run;
%ods_on();

data lsm_policy_&rdrate.;
  set lsm_mia;

  lsm_policy_nimp  = nimpute;
  lsm_policy_est   = estimate;
  lsm_policy_var   = stderr*2;
  lsm_policy_se    = stderr;
  lsm_policy_lower = ifn( lclmean ne ., lclmean, estimate - probit(0.975)*stderr);  *** IF NO MISSING DATA LCLMEAN IS MISSING ***;
  lsm_policy_upper = ifn( uclmean ne ., uclmean, estimate + probit(0.975)*stderr);  *** IF NO MISSING DATA UCLMEAN IS MISSING ***;

  if fitn=1 then do;
    select(visitn);
     when (1) miseed = &seed1.;
     when (2) miseed = &seed2.;
     when (3) miseed = &seed3.;
     when (4) miseed = &seed4.;
     when (5) miseed = &seed5.;
     otherwise;
    end;
  end;
  else do;
    miseed = &seed.;
  end;
  keep rdrate sim_run visitn groupn fit: lsm_: miseed;

run;


********************************************************************;
*** COMBINE DIF RESULTS USING RUBINS RULES                       ***;
********************************************************************;

proc sort data = dif;
  by rdrate sim_run fitn fit groupn _groupn visitn imputation;
run;

%ods_off();
proc mianalyze data = dif;
  by rdrate sim_run fitn fit groupn _groupn visitn;
  modeleffects estimate;
  stderr stderr;
  ods output parameterestimates = dif_mia;
run;
%ods_on();

data dif_policy_&rdrate.;
  set dif_mia;

  dif_policy_est   = estimate;
  dif_policy_var   = stderr*2;
  dif_policy_se    = stderr;
  dif_policy_lower = ifn( lclmean ne ., lclmean, estimate - probit(0.975)*stderr);  *** IF NO MISSING DATA LCLMEAN IS MISSING ***;
  dif_policy_upper = ifn( uclmean ne ., uclmean, estimate + probit(0.975)*stderr);  *** IF NO MISSING DATA UCLMEAN IS MISSING ***;
  _visitn           = visitn;

  if fitn=2 then do;
    select(visitn);
     when (1) miseed = &seed1.;
     when (2) miseed = &seed2.;
     when (3) miseed = &seed3.;
     when (4) miseed = &seed4.;
     when (5) miseed = &seed5.;
     otherwise;
    end;
  end;
  else do;
    miseed = &seed.;
  end;
  keep rdrate sim_run groupn _groupn visitn _visitn fit: dif_: miseed;

run;

%mend run_mi;


********************************************************************;
*** CALL ONCE FOR EACH RDRATE                                    ***;
********************************************************************;

%run_mi(rdrate=1, seed=1111, nimp=25);
%run_mi(rdrate=2, seed=2222, nimp=25);
%run_mi(rdrate=3, seed=3333, nimp=25);
%run_mi(rdrate=4, seed=4444, nimp=25);
%run_mi(rdrate=5, seed=5555, nimp=25);
%run_mi(rdrate=6, seed=6666, nimp=25);


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

proc sort data = lsm_policy (keep = rdrate sim_run groupn visitn fitn fit lsm_policy_: miseed)
          out  = results.&scenario._mi2_lsm_part&part.;
  by rdrate sim_run groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn fitn fit dif_policy_: miseed)
          out  = results.&scenario._mi2_dif_part&part.;
  by rdrate sim_run visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY AND REMOVE TEMP MI DATA FROM RESULTS    ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: mi: lsm: dif: ;
quit;





















/*******************************************************************
| Name       : dnar_drop_cr_part1.sas
| Purpose    : Code to fit copy reference to sim study data 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 01MAY22 
|-------------------------------------------------------------------
| Notes: CR
|
| Standard copy reference programmed from scratch using the 
| same basic process as the LSHTM/GSK "five macros".
|
| Basic Idea:
| 1. Parameter estimation model fitted to all the observed data
|    using a basic MAR assuption. This model is used to get draws
|    from the posterior of this model.
| 2. Imputation model created using the draws from the posterior
|    combined with the subjects covariates and values.   
| 3. Analysis model run on each completed imputation dataset and
|    then combined via rubins rules.
|
| SAS Code:
| 1. The parameter estimation is fitted and the draws from the 
|    posterior are obtained using BGLIMM. Due to MCMC optimization
|    this is much quicker than using PROC MCMC.
| 2. The creation of the copy reference means and variances 
|    via the imputation model is done using ARRAYS and these
|    are fed into an FCMP function COND_IMPUTE_MVN to do the MI.
| 3. The analysis model is a simple ANCOVA analysis by VISITN.
|
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%let scenario = dnar_drop;
%let part     = 1;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/rbi_tools.sas";


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


%macro run_mi(rdrate=, seed=, nimp=);


********************************************************************;
*** FIT PARAMETER ESTIMATION MODEL TO ON TREATMENT DATA ONLY     ***;
********************************************************************;

%ods_off();
proc bglimm data = ds_long (where = (rdrate = &rdrate. and disc = 0))
            outpost = ps 
            seed = &seed.
            nbi = 1000
            nmc = &nimp.
            thin = 1
            noclprint;
  by rdrate sim_run;
  class groupn visitn subjid;
  model response = visitn*groupn visitn*baseline_var / noint coeffprior=normal;
  repeated visitn / subject = subjid type = un group = groupn;
run;
%ods_on();


********************************************************************;
*** RENAME PARAMETER VARIABLES TO MAKE IMPUTATION MODEL EASIER   ***;
********************************************************************;

proc transpose data = ps
               out  = ps_long (rename=(col1 = value));
  by rdrate sim_run iteration;
  var groupn_: baseline_var_: residual_un_:;
run;

data ps_long;
  set ps_long;
  by rdrate sim_run;
  
  imputation = iteration - 1000;
  
  if index(lowcase(_name_), "groupn_") > 0 and index(lowcase(_name_), "visitn_") > 0 then do;
     varname = "tt"||scan(_name_, 2, "_")||scan(_name_, 4, "_");
  end;
  else if index(lowcase(_name_), "baseline_var_") > 0 and index(lowcase(_name_), "visitn_") > 0 then do;
    varname = "bt"||scan(_name_, 4, "_");
  end;
  else if index(lowcase(_name_), "residual_un_") > 0 then do;
    varname = "vc"||scan(_name_, 6, "_")||"_"||scan(_name_, 3, "_")||scan(_name_, 4, "_");  
  end;

run;

proc transpose data = ps_long
               out  = draws;
  by rdrate sim_run imputation;
  var value;
  id varname;
run;


********************************************************************;
*** SPLIT OUT COMPLETE AND MISSING DATA AND MERGE ON DRAWS       ***;
********************************************************************;

data ds1_comp
     ds1_miss;
  set ds_wide (where = (rdrate = &rdrate.));
  if withtime = 5 then output ds1_comp;
  else output ds1_miss;
run;

data ds1_comp;
  set ds1_comp;
  do imputation = 1 to &nimp.;
    output;
  end;
run;

proc sql noprint;
  create table ds1_imps (drop = b_:) as
    select *
    from ds1_miss, draws (rename = (rdrate = b_rdrate sim_run = b_sim_run)) 
    where rdrate = b_rdrate and sim_run = b_sim_run 
    order by rdrate, sim_run, groupn, subjid, imputation;
quit;


********************************************************************;
*** CREATE JUMP TO REFERENCE MU AND VC AND CONDITIONALLY IMPUTE  ***;
********************************************************************;

options nonotes;

data ds1_cr;
  set ds1_imps;
  
  *** ARRAYS TO HOLD DATA***;
  array y[5] y1-y5;
  array d[5] d1-d5;
  array i[5] i1-i5;
  
  *** ARRAYS TO HOLD IMPUTATION MODEL PARAMETERS ***;
  array tt[2,5] tt11-tt15 tt21-tt25; 
  array bt[5] bt1-bt5;
  array vc[2,15] vc1_11 vc1_21-vc1_22 vc1_31-vc1_33 vc1_41-vc1_44 vc1_51-vc1_55
                 vc2_11 vc2_21-vc2_22 vc2_31-vc2_33 vc2_41-vc2_44 vc2_51-vc2_55;
  
  *** ARRAYS TO HOLD IMPUTATION MODEL MEAN AND VARIANCE ***;
  array m[5] m1-m5;
  array v[5,5] v11-v15 v21-v25 v31-v35 v41-v45 v51-v55;
  
  *** IMPUTATION MODEL FOR MEAN ***;
  do j = 1 to 5;
    m[j] = tt[2,j] + bt[j]*baseline_var;  *** REF MEAN USED FOR ALL TIMEPOINTS ***;
  end;
    
  *** IMPUTATION MODEL FOR VARIANCE ***;  
  do j = 1 to 5;
  do k = 1 to j;
    v[j,k] = vc[2,(j*(j-1)/2)+k];  *** FILL LOWER V FROM REF VC ***;
    v[k,j] = v[j,k];               *** CREATE SYMMETRIC V MAT ***;
  end;
  end;
  
  *** USE MEAN AND V FOR CONDITIONAL IMPUTATION ***;
  call cond_impute_mvn(y, m, v, i);
      
  drop tt: bt: vc: m1-m5 v11-v15 v21-v25 v31-v35 v41-v45 v51-v55;   
      
run;

options notes;


********************************************************************;
*** RESTACK COMPLETE AND IMPUTED DATA TOGETHER AND SORT          ***;
********************************************************************;

data mi;
  set ds1_comp
      ds1_cr;
  by rdrate sim_run;
run;

proc sort data = mi;
  by rdrate sim_run imputation groupn subjid;
run;


********************************************************************;
*** TRANSPOSE MI DATA BACK INTO LONG FORMAT FOR MIXED            ***;
********************************************************************;

data mi_long;
  set mi;
  by rdrate sim_run imputation groupn subjid;

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
      lsm_policy_6
      ;
run;

data dif_policy;
  set dif_policy_1
      dif_policy_2
      dif_policy_3
      dif_policy_4
      dif_policy_5
      dif_policy_6
      ;
run;

********************************************************************;
*** ORGANISE RESULTS AND OUTPUT TO PERM LOCATION                 ***;
********************************************************************;

proc sort data = lsm_policy (keep = rdrate sim_run groupn visitn lsm_policy_: miseed)
          out  = results.&scenario._cr_lsm_part&part.;
  by rdrate sim_run groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn dif_policy_: miseed)
          out  = results.&scenario._cr_dif_part&part.;
  by rdrate sim_run visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY AND REMOVE TEMP MI DATA FROM RESULTS    ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: ps: draws mi: lsm: dif: functions ;
quit;



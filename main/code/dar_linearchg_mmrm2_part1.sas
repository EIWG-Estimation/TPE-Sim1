/*******************************************************************
| Name       : dar_linearchg_mmrm2_part1.sas
| Purpose    : Code to fit MMRM with TDC and OM to sim study data 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 15SEP21 
|-------------------------------------------------------------------
| Notes: MMRM2
|
| This is an MMRM fitted to CHANGE with time dependent covariate 
| for On/Off treatment with manually weighted observed margins. 
|     
| The use of the time dependent covariate splits the estimation 
| into On and Off treatment states.
|
| Treatment Policy Estimates are then created as a weighted average
| over the proportions On and Off at each time.
|
| The variability is then corrected for the fact that the
| proportions are estimates of the true proportion. 
|
| Variance correction assumes that proportions on/off are not
| dependent on the expected values. 
|
| If the MMRM2 model is non-estimable for any of the treatment by 
| time by discontinuation status levels then these results are not
| used and the simpler naive MMRM1 model is fitted to that simulated
| trial. This is a pre-specified and principled analysis.
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

proc sort data = ds_long;
    by rdrate sim_run patmaxn patmax visitn groupn subjid;
run;


********************************************************************;
*** FIT MMRM2 USING DISCCODE COVARIATE AND CORRECT BASELINE MEAN ***;
********************************************************************;

%macro mmrm2(rdrate=);

  %ods_off(notes=N);
   
  %let nsims = 0;    *** IN CASE THERE ARE NO CASES IN DATA ***;
   
  *** COUNT SIMS ***;
  data ds_mmrm2;
    set ds_long (where = (rdrate = &rdrate.)) end = eof;
    retain sim_count 0;
    by rdrate sim_run patmaxn patmax;
    if first.sim_run then sim_count + 1;
	if eof then call symput("nsims", put(sim_count,best.));
  run;

  %do ii = 1 %to &nsims.;

  *** SQL TO CREATE MACRO VARIABLE FOR BASELINE MEAN ***;
  proc sql noprint;
    select mean(baseline_var) into :mbl
    from ds_mmrm2
    where sim_count = &ii. and visitn = 1;
  quit;

  *** MMRM2 WITH LSM AT CORRECT BASELINE VALUE ***;
  proc mixed data = ds_mmrm2 (where = (sim_count = &ii.));
    by rdrate sim_run patmaxn patmax;
    class visitn groupn disccode subjid;
    model change = visitn*groupn*disccode visitn*baseline_var / noint ddfm=kr;
    repeated visitn / subject=subjid type=un;
    lsmeans visitn*groupn*disccode / at baseline_var = &mbl. cov;
    ods output lsmeans = lsm_sc&ii.;
  run;

  *** STACK LSM RESULTS FOR EACH SIM_RUN ***; 
  data lsm_r&rdrate.;
    set %if &ii. ne 1 %then %do; lsm_r&rdrate. %end;
        lsm_sc&ii.;
    by rdrate sim_run patmaxn patmax visitn groupn disccode;
  run;

  %end;
  %ods_on();

%mend;

%mmrm2(rdrate=1);
%mmrm2(rdrate=2);
%mmrm2(rdrate=3);
%mmrm2(rdrate=4);
%mmrm2(rdrate=5);
%mmrm2(rdrate=6);


********************************************************************;
*** STACK RESULTS AND CLEAN UP WORK AREA                         ***;
********************************************************************;

data lsm_mmrm2;
  set lsm_r:;
  by rdrate sim_run patmaxn patmax visitn groupn disccode;
run; 

%ods_off();
proc datasets lib = work nolist;
  delete lsm_sc: lsm_r: ds_mmrm2;
quit;
%ods_on();


********************************************************************;
*** CALCULATE PROPORTIONS BY DISCCODE AND TRANSPOSE              ***;
********************************************************************;

proc freq data = ds_long noprint;
  by rdrate sim_run patmaxn patmax visitn groupn;
  tables disccode / out = mmrm2_om;
run;

data lsm_mmrm2_om;
  set mmrm2_om;
  by rdrate sim_run patmaxn patmax visitn groupn disccode;
  retain row pat 0 n_pat1-n_pat2 .;

  array n_pat[2] n_pat1-n_pat2;

  if first.groupn then do;
    pat = 0;
    do i = 1 to 2;
      n_pat[i] = .;
    end;
  end;

  if first.disccode then pat = pat + 1;
  n_pat[pat] = count;                     *** CREATE THE TRANSPOSED COUNTS ***;

  if last.groupn then output;

  keep rdrate sim_run patmaxn patmax visitn groupn n_pat:;
 
run;


********************************************************************;
*** TRANSPOSE THE LSM AND VC ESTIMATES                           ***;
********************************************************************;

data lsm_mmrm2_wide;
  set lsm_mmrm2;
  by rdrate sim_run patmaxn patmax visitn groupn disccode;
  retain row pat 0 lsm_pat1-lsm_pat2
         lsm_vc11 
         lsm_vc21-lsm_vc22 .;

  array covs[20] cov1-cov20;

  array lsm_pat[2] lsm_pat1-lsm_pat2;
  array lsm_vc[3] lsm_vc11  lsm_vc21-lsm_vc22; 

  if first.sim_run then row = 0;
  if first.groupn then do;
    pat = 0;
    do i = 1 to 2;
      lsm_pat[i] = .;
    end;
    do i = 1 to 3;
      lsm_vc[i] = .;
    end;
  end;

  if first.disccode then pat = pat + 1;
  row = row + 1;

  *** CREATE THE TRANSPOSED LSM ***;
  lsm_pat[pat] = estimate;  


  *** PICK OUT THE PATTERN VARIANCES AND COVARIANCE ***;
  select (pat);
    when (1) do; lsm_vc11 = covs[row-0]; end;
    when (2) do; lsm_vc21 = covs[row-1]; lsm_vc22 = covs[row-0]; end;
    otherwise;
  end;
    
  if last.groupn then output;
  keep rdrate sim_run patmaxn patmax visitn groupn lsm_pat: lsm_vc:;
 
run;


********************************************************************;
*** MERGE PROPORTIONS LSMS AND VCS AND CALCULATE POLICY          ***;
********************************************************************;

data lsm_policy;
  merge lsm_mmrm2_om
        lsm_mmrm2_wide;
  by rdrate sim_run patmaxn patmax visitn groupn;
  
  ntotal = sum(of n_pat1-n_pat2);
  npats  = 2-nmiss(lsm_pat1,lsm_pat2);

  if npats = 1 then do;         *** FORMULA FOR POLICY AND VARIANCE ADJUSTMENT BASED ON PATTERNS ***;
  
    lsm_policy_est     = lsm_pat1;  *** IF ONLY ONE DISCCODE PAT POSSIBLE NO WEIGHING NEEDED ***;
    lsm_policy_var     = lsm_vc11;  
    lsm_policy_var_adj = lsm_vc11;  *** IF ONLY ONE DISCCODE PAT POSSIBLE NO DELTA METHOD ADJUSTMENT NEEDED ***;
    
  end;
  else if npats = 2 then do;
    
    p_pat1 = n_pat1 / ntotal;
    p_pat2 = n_pat2 / ntotal;  

    lsm_policy_est = p_pat1*lsm_pat1 + p_pat2*lsm_pat2;

    vc_p11  = (p_pat1 * (1-p_pat1))  / ntotal;  *** VAR OF PROPORTION FOR PAT1 ***;    

    delta12 = lsm_pat1 - lsm_pat2;    *** DIF BETWEEN LSM PAT1 AND PAT2 ***;

    lsm_policy_var = 1*lsm_vc11*p_pat1*p_pat1 + 
                     2*lsm_vc21*p_pat2*p_pat1 + 1*lsm_vc22*p_pat2*p_pat2 ;

    lsm_policy_var_adj = lsm_policy_var + 
                         1*vc_p11*delta12*delta12;

  end;

  lsm_policy_se    = sqrt(lsm_policy_var); 
  lsm_policy_lower = lsm_policy_est - probit(0.975)*lsm_policy_se;
  lsm_policy_upper = lsm_policy_est + probit(0.975)*lsm_policy_se;

  lsm_policy_se_adj    = sqrt(lsm_policy_var_adj); 
  lsm_policy_lower_adj = lsm_policy_est - probit(0.975)*lsm_policy_se_adj;
  lsm_policy_upper_adj = lsm_policy_est + probit(0.975)*lsm_policy_se_adj;

  length fit $50;
  if patmaxn in (2 3 4 5 6) then do;
    fitn = 1;
	fit  = "2-pattern MMRM (MMRM2)";
  end;
  else if patmaxn = 1 then do;
    fitn = 2;
	fit  = "1-pattern MMRM (MMRM1)";
  end;

  keep rdrate sim_run patmaxn patmax fitn fit visitn groupn ntotal lsm_policy_:;
  
run;


********************************************************************;
*** CALC DIFS IN LSM                                             ***;
********************************************************************;

data dif_policy;
  merge lsm_policy (where=(groupn=1) rename=(ntotal = n1 lsm_policy_est=lsm_policy1 lsm_policy_var=lsm_policy_var1 lsm_policy_var_adj=lsm_policy_var_adj1))
        lsm_policy (where=(groupn=2) rename=(ntotal = n2 lsm_policy_est=lsm_policy2 lsm_policy_var=lsm_policy_var2 lsm_policy_var_adj=lsm_policy_var_adj2));
  by rdrate sim_run patmaxn patmax visitn;

  groupn  = 1;
  _groupn = 2;
  _visitn = visitn;

  dif_policy_est   = lsm_policy1 - lsm_policy2;
  
  dif_policy_var   = lsm_policy_var1 + lsm_policy_var2;
  dif_policy_se    = sqrt( dif_policy_var );
  dif_policy_lower = dif_policy_est - probit(0.975)*dif_policy_se;
  dif_policy_upper = dif_policy_est + probit(0.975)*dif_policy_se;

  dif_policy_var_adj   = lsm_policy_var_adj1 + lsm_policy_var_adj2;
  dif_policy_se_adj    = sqrt( dif_policy_var_adj );
  dif_policy_lower_adj = dif_policy_est - probit(0.975)*dif_policy_se_adj;
  dif_policy_upper_adj = dif_policy_est + probit(0.975)*dif_policy_se_adj;

  keep rdrate sim_run patmaxn patmax fitn fit visitn _visitn groupn _groupn n1 n2 dif_policy_:;

run;


********************************************************************;
*** ORGANISE RESULTS AND OUTPUT TO PERM LOCATION                 ***;
********************************************************************;

proc sort data = lsm_policy (keep = rdrate sim_run groupn visitn patmaxn patmax fitn fit lsm_policy_:)
          out  = results.&scenario._mmrm2_lsm_part&part.;
  by rdrate sim_run patmaxn patmax fitn fit groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn patmaxn patmax fitn fit dif_policy_:)
          out  = results.&scenario._mmrm2_dif_part&part.;
  by rdrate sim_run patmaxn patmax fitn fit visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY                                         ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: lsm: dif: mmrm: ;
quit;


/*******************************************************************
| Name       : dar_drop_mmrm3_part1.sas
| Purpose    : Code to fit MMRM with TDC and OM to sim study data 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 15APR22 
|-------------------------------------------------------------------
| Notes: MMRM3 Complete Step Down
|
| An MMRM fitted to CHANGE with time dependent covariate for pattern 
| of discontinuation with manually weighted observed margins. 
|     
| The use of the time dependent covariate splits the estimation 
| into treatment discontinuation patterns.
|
| Treatment Policy Estimates are then created as a weighted average
| over the proportions in each pattern at each time.
|
| The variability is then corrected for the fact that the
| proportions are estimates of the true proportion. 
|
| Variance correction assumes that proportions are not
| dependent on the expected values. 
|
| If the MMRM3 model is non-estimable for any of the treatment by 
| time by discontinuation pattern levels then these results are not
| used and a step down algorithm is applied. The step down combines
| any pattern that is non-estimable with the preceeding pattern
| unless that is the first discontinuation pattern which is
| is combined with the second discontinuation pattern.
| 
| The last pattern at any visit is always the on-treatment pattern.
| 
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%let scenario = dar_drop;
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
*** FIT MMRM3 USING PATTERN COVARIATE AND CORRECT BASELINE MEAN  ***;
********************************************************************;

%macro mmrm3(rdrate=);

  %ods_off(notes=N);
   
  %let nsims = 0;    *** IN CASE THERE ARE NO CASES IN DATA ***;
   
  *** COUNT SIMS ***;
  data ds_mmrm3;
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
    from ds_mmrm3
    where sim_count = &ii. and visitn = 1;
  quit;

  *** MMRM3 WITH LSM AT CORRECT BASELINE VALUE ***;
  proc mixed data = ds_mmrm3 (where = (sim_count = &ii.));
    by rdrate sim_run patmaxn patmax;
    class visitn groupn patcode subjid;
    model change = visitn*groupn*patcode visitn*baseline_var / noint ddfm=kr;
    repeated visitn / subject=subjid type=un;
    lsmeans visitn*groupn*patcode / at baseline_var = &mbl. cov;
    ods output lsmeans = lsm_sc&ii.;
  run;

  *** STACK LSM RESULTS FOR EACH SIM_RUN ***; 
  data lsm_r&rdrate.;
    set %if &ii. ne 1 %then %do; lsm_r&rdrate. %end;
        lsm_sc&ii.;
    by rdrate sim_run patmaxn patmax visitn groupn patcode;
  run;

  %end;
  %ods_on();

%mend;

%mmrm3(rdrate=1);
%mmrm3(rdrate=2);
%mmrm3(rdrate=3);
%mmrm3(rdrate=4);
%mmrm3(rdrate=5);
%mmrm3(rdrate=6);


********************************************************************;
*** STACK RESULTS AND CLEAN UP WORK AREA                         ***;
********************************************************************;

data lsm_mmrm3;
  set lsm_r:;
  by rdrate sim_run patmaxn patmax visitn groupn patcode;
run; 

%ods_off();
proc datasets lib = work nolist;
  delete lsm_sc: lsm_r: ds_mmrm3;
quit;
%ods_on();


********************************************************************;
*** CALCULATE PROPORTIONS BY PATTERN AND TRANSPOSE               ***;
********************************************************************;

proc freq data = ds_long noprint;
  by rdrate sim_run patmaxn patmax visitn groupn;
  tables patcode / out = mmrm3_om;
run;

data lsm_mmrm3_om;
  set mmrm3_om;
  by rdrate sim_run patmaxn patmax visitn groupn patcode;
  retain row pat 0 n_pat1-n_pat6 .;

  array n_pat[6] n_pat1-n_pat6;

  if first.groupn then do;
    pat = 0;
    do i = 1 to 6;
      n_pat[i] = .;
    end;
  end;

  if first.patcode then pat = pat + 1;
  n_pat[pat] = count;                     *** CREATE THE TRANSPOSED COUNTS ***;

  if last.groupn then output;

  keep rdrate sim_run patmaxn patmax visitn groupn n_pat:;
 
run;


********************************************************************;
*** TRANSPOSE THE LSM AND VC ESTIMATES                           ***;
********************************************************************;

data lsm_mmrm3_wide;
  set lsm_mmrm3;
  by rdrate sim_run patmaxn patmax visitn groupn patcode;
  retain row pat 0 lsm_pat1-lsm_pat6
         lsm_vc11 
         lsm_vc21-lsm_vc22 
         lsm_vc31-lsm_vc33
         lsm_vc41-lsm_vc44
         lsm_vc51-lsm_vc55
         lsm_vc61-lsm_vc66 .;

  array covs[40] cov1-cov40;

  array lsm_pat[6] lsm_pat1-lsm_pat6;
  array lsm_vc[21] lsm_vc11  lsm_vc21-lsm_vc22  lsm_vc31-lsm_vc33  lsm_vc41-lsm_vc44  lsm_vc51-lsm_vc55  lsm_vc61-lsm_vc66; 

  if first.sim_run then row = 0;
  if first.groupn then do;
    pat = 0;
    do i = 1 to 6;
      lsm_pat[i] = .;
    end;
    do i = 1 to 21;
      lsm_vc[i] = .;
    end;
  end;

  if first.patcode then pat = pat + 1;
  row = row + 1;

  *** CREATE THE TRANSPOSED LSM ***;
  lsm_pat[pat] = estimate;  


  *** PICK OUT THE PATTERN VARIANCES AND COVARIANCE ***;
  select (pat);
    when (1) do; lsm_vc11 = covs[row-0]; end;
    when (2) do; lsm_vc21 = covs[row-1]; lsm_vc22 = covs[row-0]; end;
    when (3) do; lsm_vc31 = covs[row-2]; lsm_vc32 = covs[row-1]; lsm_vc33 = covs[row-0]; end;
    when (4) do; lsm_vc41 = covs[row-3]; lsm_vc42 = covs[row-2]; lsm_vc43 = covs[row-1]; lsm_vc44 = covs[row-0];end;
    when (5) do; lsm_vc51 = covs[row-4]; lsm_vc52 = covs[row-3]; lsm_vc53 = covs[row-2]; lsm_vc54 = covs[row-1]; lsm_vc55 = covs[row-0];end;
    when (6) do; lsm_vc61 = covs[row-5]; lsm_vc62 = covs[row-4]; lsm_vc63 = covs[row-3]; lsm_vc64 = covs[row-2]; lsm_vc65 = covs[row-1]; lsm_vc66 = covs[row-0];end;
    otherwise;
  end;
    
  if last.groupn then output;
  keep rdrate sim_run patmaxn patmax visitn groupn lsm_pat: lsm_vc:;
 
run;


********************************************************************;
*** MERGE PROPORTIONS LSMS AND VCS AND CALCULATE POLICY          ***;
********************************************************************;

data lsm_policy;
  merge lsm_mmrm3_om
        lsm_mmrm3_wide;
  by rdrate sim_run patmaxn patmax visitn groupn;
  
  ntotal = sum(of n_pat1-n_pat6);
  npats  = 6-nmiss(lsm_pat1,lsm_pat2,lsm_pat3,lsm_pat4,lsm_pat5,lsm_pat6);

  if npats = 1 then do;         *** FORMULA FOR POLICY AND VARIANCE ADJUSTMENT BASED ON PATTERNS ***;
  
    lsm_policy_est     = lsm_pat1;  *** IF ONLY ONE PAT POSSIBLE NO WEIGHING NEEDED ***;
    lsm_policy_var     = lsm_vc11;  
    lsm_policy_var_adj = lsm_vc11;  *** IF ONLY ONE PAT POSSIBLE NO DELTA METHOD ADJUSTMENT NEEDED ***;
    
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
  else if npats = 3 then do;
    
    p_pat1 = n_pat1 / ntotal;
    p_pat2 = n_pat2 / ntotal;
    p_pat3 = n_pat3 / ntotal;  

    lsm_policy_est = p_pat1*lsm_pat1 + p_pat2*lsm_pat2 + p_pat3*lsm_pat3;

    vc_p11  = (p_pat1 * (1-p_pat1))  / ntotal;  *** VAR OF PROPORTION FOR PAT1 ***;
    vc_p22  = (p_pat2 * (1-p_pat2))  / ntotal;  *** VAR OF PROPORTION FOR PAT2 ***;

    vc_p12 = (-1 * p_pat1 * p_pat2) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT2 ***;
      
    delta13 = lsm_pat1 - lsm_pat3;    *** DIF BETWEEN LSM PAT1 AND PAT3 ***;
    delta23 = lsm_pat2 - lsm_pat3;    *** DIF BETWEEN LSM PAT2 AND PAT3 ***;

    lsm_policy_var = 1*lsm_vc11*p_pat1*p_pat1 + 
                     2*lsm_vc21*p_pat2*p_pat1 + 1*lsm_vc22*p_pat2*p_pat2 +
                     2*lsm_vc31*p_pat3*p_pat1 + 2*lsm_vc32*p_pat3*p_pat2 + 1*lsm_vc33*p_pat3*p_pat3;

    lsm_policy_var_adj = lsm_policy_var + 
                         1*vc_p11*delta13*delta13 + 
                         2*vc_p12*delta13*delta23 + 1*vc_p22*delta23*delta23;

  end;
  else if npats = 4 then do;
  
    p_pat1 = n_pat1 / ntotal;
    p_pat2 = n_pat2 / ntotal;
    p_pat3 = n_pat3 / ntotal;
    p_pat4 = n_pat4 / ntotal;  

    lsm_policy_est = p_pat1*lsm_pat1 + p_pat2*lsm_pat2 + p_pat3*lsm_pat3 + p_pat4*lsm_pat4;

    vc_p11  = (p_pat1 * (1-p_pat1))  / ntotal;  *** VAR OF PROPORTION FOR PAT1 ***;
    vc_p22  = (p_pat2 * (1-p_pat2))  / ntotal;  *** VAR OF PROPORTION FOR PAT2 ***;
    vc_p33  = (p_pat3 * (1-p_pat3))  / ntotal;  *** VAR OF PROPORTION FOR PAT3 ***;

    vc_p12 = (-1 * p_pat1 * p_pat2) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT2 ***;
    vc_p13 = (-1 * p_pat1 * p_pat3) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT3 ***;
    vc_p23 = (-1 * p_pat2 * p_pat3) / ntotal;   *** COV OF PROPORTION FOR PAT2 AND PAT3 ***;
      
    delta14 = lsm_pat1 - lsm_pat4;    *** DIF BETWEEN LSM PAT1 AND PAT4 ***;
    delta24 = lsm_pat2 - lsm_pat4;    *** DIF BETWEEN LSM PAT2 AND PAT4 ***;
    delta34 = lsm_pat3 - lsm_pat4;    *** DIF BETWEEN LSM PAT3 AND PAT4 ***;

    lsm_policy_var = 1*lsm_vc11*p_pat1*p_pat1 + 
                     2*lsm_vc21*p_pat2*p_pat1 + 1*lsm_vc22*p_pat2*p_pat2 +
                     2*lsm_vc31*p_pat3*p_pat1 + 2*lsm_vc32*p_pat3*p_pat2 + 1*lsm_vc33*p_pat3*p_pat3 +
                     2*lsm_vc41*p_pat4*p_pat1 + 2*lsm_vc42*p_pat4*p_pat2 + 2*lsm_vc43*p_pat4*p_pat3 + 1*lsm_vc44*p_pat4*p_pat4;

    lsm_policy_var_adj = lsm_policy_var + 
                         1*vc_p11*delta14*delta14 + 
                         2*vc_p12*delta14*delta24 + 1*vc_p22*delta24*delta24 +
                         2*vc_p13*delta14*delta34 + 2*vc_p23*delta24*delta34 + 1*vc_p33*delta34*delta34; 

  end;
  else if npats = 5 then do;

    p_pat1 = n_pat1 / ntotal;
    p_pat2 = n_pat2 / ntotal;
    p_pat3 = n_pat3 / ntotal;
    p_pat4 = n_pat4 / ntotal;
    p_pat5 = n_pat5 / ntotal;  

    lsm_policy_est = p_pat1*lsm_pat1 + p_pat2*lsm_pat2 + p_pat3*lsm_pat3 + p_pat4*lsm_pat4 + p_pat5*lsm_pat5;

    vc_p11  = (p_pat1 * (1-p_pat1))  / ntotal;  *** VAR OF PROPORTION FOR PAT1 ***;
    vc_p22  = (p_pat2 * (1-p_pat2))  / ntotal;  *** VAR OF PROPORTION FOR PAT2 ***;
    vc_p33  = (p_pat3 * (1-p_pat3))  / ntotal;  *** VAR OF PROPORTION FOR PAT3 ***;
    vc_p44  = (p_pat4 * (1-p_pat4))  / ntotal;  *** VAR OF PROPORTION FOR PAT4 ***;

    vc_p12 = (-1 * p_pat1 * p_pat2) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT2 ***;
    vc_p13 = (-1 * p_pat1 * p_pat3) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT3 ***;
    vc_p23 = (-1 * p_pat2 * p_pat3) / ntotal;   *** COV OF PROPORTION FOR PAT2 AND PAT3 ***;
    vc_p14 = (-1 * p_pat1 * p_pat4) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT4 ***;
    vc_p24 = (-1 * p_pat2 * p_pat4) / ntotal;   *** COV OF PROPORTION FOR PAT2 AND PAT4 ***;
    vc_p34 = (-1 * p_pat3 * p_pat4) / ntotal;   *** COV OF PROPORTION FOR PAT3 AND PAT4 ***;
      
    delta15 = lsm_pat1 - lsm_pat5;    *** DIF BETWEEN LSM PAT1 AND PAT5 ***;
    delta25 = lsm_pat2 - lsm_pat5;    *** DIF BETWEEN LSM PAT2 AND PAT5 ***;
    delta35 = lsm_pat3 - lsm_pat5;    *** DIF BETWEEN LSM PAT3 AND PAT5 ***;
    delta45 = lsm_pat4 - lsm_pat5;    *** DIF BETWEEN LSM PAT4 AND PAT5 ***;      

    lsm_policy_var = 1*lsm_vc11*p_pat1*p_pat1 + 
                     2*lsm_vc21*p_pat2*p_pat1 + 1*lsm_vc22*p_pat2*p_pat2 +
                     2*lsm_vc31*p_pat3*p_pat1 + 2*lsm_vc32*p_pat3*p_pat2 + 1*lsm_vc33*p_pat3*p_pat3 +
                     2*lsm_vc41*p_pat4*p_pat1 + 2*lsm_vc42*p_pat4*p_pat2 + 2*lsm_vc43*p_pat4*p_pat3 + 1*lsm_vc44*p_pat4*p_pat4 +
                     2*lsm_vc51*p_pat5*p_pat1 + 2*lsm_vc52*p_pat5*p_pat2 + 2*lsm_vc53*p_pat5*p_pat3 + 2*lsm_vc54*p_pat5*p_pat4 + 1*lsm_vc55*p_pat5*p_pat5; 
  
    lsm_policy_var_adj = lsm_policy_var + 
                         1*vc_p11*delta15*delta15 + 
                         2*vc_p12*delta15*delta25 + 1*vc_p22*delta25*delta25 +
                         2*vc_p13*delta15*delta35 + 2*vc_p23*delta25*delta35 + 1*vc_p33*delta35*delta35 + 
                         2*vc_p14*delta15*delta45 + 2*vc_p24*delta25*delta45 + 2*vc_p34*delta35*delta45 + 1*vc_p44*delta45*delta45; 

  end;
  else if npats = 6 then do;
 
    p_pat1 = n_pat1 / ntotal;
    p_pat2 = n_pat2 / ntotal;
    p_pat3 = n_pat3 / ntotal;
    p_pat4 = n_pat4 / ntotal;
    p_pat5 = n_pat5 / ntotal;  
    p_pat6 = n_pat6 / ntotal;  

    lsm_policy_est = p_pat1*lsm_pat1 + p_pat2*lsm_pat2 + p_pat3*lsm_pat3 + p_pat4*lsm_pat4 + p_pat5*lsm_pat5 + p_pat6*lsm_pat6;

    vc_p11  = (p_pat1 * (1-p_pat1))  / ntotal;  *** VAR OF PROPORTION FOR PAT1 ***;
    vc_p22  = (p_pat2 * (1-p_pat2))  / ntotal;  *** VAR OF PROPORTION FOR PAT2 ***;
    vc_p33  = (p_pat3 * (1-p_pat3))  / ntotal;  *** VAR OF PROPORTION FOR PAT3 ***;
    vc_p44  = (p_pat4 * (1-p_pat4))  / ntotal;  *** VAR OF PROPORTION FOR PAT4 ***;
    vc_p55  = (p_pat5 * (1-p_pat5))  / ntotal;  *** VAR OF PROPORTION FOR PAT5 ***;

    vc_p12 = (-1 * p_pat1 * p_pat2) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT2 ***;
    vc_p13 = (-1 * p_pat1 * p_pat3) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT3 ***;
    vc_p23 = (-1 * p_pat2 * p_pat3) / ntotal;   *** COV OF PROPORTION FOR PAT2 AND PAT3 ***;
    vc_p14 = (-1 * p_pat1 * p_pat4) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT4 ***;
    vc_p24 = (-1 * p_pat2 * p_pat4) / ntotal;   *** COV OF PROPORTION FOR PAT2 AND PAT4 ***;
    vc_p34 = (-1 * p_pat3 * p_pat4) / ntotal;   *** COV OF PROPORTION FOR PAT3 AND PAT4 ***;
    vc_p15 = (-1 * p_pat1 * p_pat5) / ntotal;   *** COV OF PROPORTION FOR PAT1 AND PAT5 ***;
    vc_p25 = (-1 * p_pat2 * p_pat5) / ntotal;   *** COV OF PROPORTION FOR PAT2 AND PAT5 ***;
    vc_p35 = (-1 * p_pat3 * p_pat5) / ntotal;   *** COV OF PROPORTION FOR PAT3 AND PAT5 ***;
    vc_p45 = (-1 * p_pat4 * p_pat5) / ntotal;   *** COV OF PROPORTION FOR PAT4 AND PAT5 ***;
      
    delta16 = lsm_pat1 - lsm_pat6;    *** DIF BETWEEN LSM PAT1 AND PAT6 ***;
    delta26 = lsm_pat2 - lsm_pat6;    *** DIF BETWEEN LSM PAT2 AND PAT6 ***;
    delta36 = lsm_pat3 - lsm_pat6;    *** DIF BETWEEN LSM PAT3 AND PAT6 ***;
    delta46 = lsm_pat4 - lsm_pat6;    *** DIF BETWEEN LSM PAT4 AND PAT6 ***;      
    delta56 = lsm_pat5 - lsm_pat6;    *** DIF BETWEEN LSM PAT5 AND PAT6 ***;      

    lsm_policy_var = 1*lsm_vc11*p_pat1*p_pat1 + 
                     2*lsm_vc21*p_pat2*p_pat1 + 1*lsm_vc22*p_pat2*p_pat2 +
                     2*lsm_vc31*p_pat3*p_pat1 + 2*lsm_vc32*p_pat3*p_pat2 + 1*lsm_vc33*p_pat3*p_pat3 +
                     2*lsm_vc41*p_pat4*p_pat1 + 2*lsm_vc42*p_pat4*p_pat2 + 2*lsm_vc43*p_pat4*p_pat3 + 1*lsm_vc44*p_pat4*p_pat4 +
                     2*lsm_vc51*p_pat5*p_pat1 + 2*lsm_vc52*p_pat5*p_pat2 + 2*lsm_vc53*p_pat5*p_pat3 + 2*lsm_vc54*p_pat5*p_pat4 + 1*lsm_vc55*p_pat5*p_pat5 + 
                     2*lsm_vc61*p_pat6*p_pat1 + 2*lsm_vc62*p_pat6*p_pat2 + 2*lsm_vc63*p_pat6*p_pat3 + 2*lsm_vc64*p_pat6*p_pat4 + 2*lsm_vc65*p_pat6*p_pat5 + 1*lsm_vc66*p_pat6*p_pat6; 
  
    lsm_policy_var_adj = lsm_policy_var +
                         1*vc_p11*delta16*delta16 + 
                         2*vc_p12*delta16*delta26 + 1*vc_p22*delta26*delta26 +
                         2*vc_p13*delta16*delta36 + 2*vc_p23*delta26*delta36 + 1*vc_p33*delta36*delta36 + 
                         2*vc_p14*delta16*delta46 + 2*vc_p24*delta26*delta46 + 2*vc_p34*delta36*delta46 + 1*vc_p44*delta46*delta46 + 
                         2*vc_p15*delta16*delta56 + 2*vc_p25*delta26*delta56 + 2*vc_p35*delta36*delta56 + 2*vc_p45*delta46*delta56 + 1*vc_p55*delta56*delta56; 
  
  end;

  lsm_policy_se    = sqrt(lsm_policy_var); 
  lsm_policy_lower = lsm_policy_est - probit(0.975)*lsm_policy_se;
  lsm_policy_upper = lsm_policy_est + probit(0.975)*lsm_policy_se;

  lsm_policy_se_adj    = sqrt(lsm_policy_var_adj); 
  lsm_policy_lower_adj = lsm_policy_est - probit(0.975)*lsm_policy_se_adj;
  lsm_policy_upper_adj = lsm_policy_est + probit(0.975)*lsm_policy_se_adj;

  fitn = 6 - patmaxn + 1;
  length fit $50;
  select(fitn);
    when (1) fit = "6-pattern MMRM (MMRM3)";
    when (2) fit = "5-pattern MMRM";
    when (3) fit = "4-pattern MMRM";
    when (4) fit = "3-pattern MMRM";
    when (5) fit = "2-pattern MMRM (MMRM2)";
    when (6) fit = "1-pattern MMRM (MMRM1)";
    otherwise;
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

proc sort data = lsm_policy (keep = rdrate sim_run patmaxn patmax groupn visitn fitn fit lsm_policy_:)
          out  = results.&scenario._mmrm3_lsm_part&part.;
  by rdrate sim_run patmaxn patmax fitn fit groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run patmaxn patmax groupn _groupn visitn _visitn fitn fit dif_policy_:)
          out  = results.&scenario._mmrm3_dif_part&part.;
  by rdrate sim_run patmaxn patmax fitn fit visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY                                         ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: lsm: dif: om: mmrm: ;
quit;

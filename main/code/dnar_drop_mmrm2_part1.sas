/*******************************************************************
| Name       : dnar_drop_mmrm2_part1.sas
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

%let scenario = dnar_drop;
%let part     = 1;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";

********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

proc sql noprint;
  create table ds as
    select *, 1-ontrt as disc 
    from data.&scenario.
    where 0 lt sim_run le 5000
    order by rdrate, sim_run, groupn, subjid, visitn, disc;
quit;


********************************************************************;
*** TRANSPOSE TO WIDE FORMAT FOR MODEL FEASIBILITY CHECKING      ***;
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
*** CHECK MODEL HAS DATA FOR EACH GROUPN AND DISC IN EACH SIM    ***;
********************************************************************;

data ds2;
  set ds1;
  by rdrate sim_run groupn;
  retain row i1-i5 0 mmrm2 1;  

  array d[5] d1-d5;
  array w[5] w1-w5;
  array i[5] i1-i5;                    *** RETAINED IMPUTATION FLAGS ***;

  if first.sim_run then mmrm2 = 1;      *** RESET OVERALL MODEL FLAG AS YES ***;
  if first.groupn then do j = 1 to 5;  *** RESET BY GROUP MI CHECK FLAGS ***;
    i[j] = 0;
  end;

  do j = 1 to 5;
    if d[j]=1 and w[j]=0 then i[j] = 1; *** IF ONE SUBJECT DISCONTINUED AND DID NOT WITHDRAW IMP IS POSSIBLE FOR THAT VISIT ***;
  end;

  if last.groupn and sum(of i1-i5) lt 5 then mmrm2 = 0;  *** CHECK EACH TREATMENT IF ANY FAILS MODEL IS NOT POSSIBLE ***;
  if last.sim_run then output;                           *** OUTPUT THE MODEL FLAG FOR EACH SIMULATION ***;
  keep rdrate sim_run mmrm2;

run;


********************************************************************;
*** MERGE ON MODEL FITTING FLAG                                  ***;
********************************************************************;

data ds;
  merge ds
        ds2;
  by rdrate sim_run;
run;

proc datasets lib = work nolist;
  delete ds1 ds2;
quit;


********************************************************************;
*** FIT MMRM2 USING DISC COVARIATE AND CORRECT BASELINE MEAN     ***;
********************************************************************;

%macro mmrm2(rdrate=);

  %ods_off();

  *** SUBSET TO MMRM2 DATA AND COUNT SIMS ***;
  data ds_mmrm2;
    set ds (where = (rdrate = &rdrate. and mmrm2=1)) end = eof;
    retain sim_count 0;
    by rdrate sim_run;
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
    by rdrate sim_run;
    class visitn groupn disc subjid;
    model change = visitn*groupn*disc visitn*baseline_var / noint ddfm=kr;
    repeated visitn / subject=subjid type=un;
    lsmeans visitn*groupn*disc / at baseline_var = &mbl. cov;
    ods output lsmeans = lsm_sc&ii.;
  run;

  *** STACK LSM RESULTS FOR EACH SIM_RUN ***; 
  data lsm_r&rdrate.;
    set %if &ii. ne 1 %then %do; lsm_r&rdrate. %end;
        lsm_sc&ii.;
    by rdrate sim_run visitn groupn disc;
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
  set lsm_r1
      lsm_r2
      lsm_r3
      lsm_r4
      lsm_r5
      lsm_r6;
  by rdrate sim_run visitn groupn;
run; 

%ods_off();
proc datasets lib = work nolist;
  delete lsm_sc: lsm_r: ds_mmrm2;
quit;
%ods_on();

********************************************************************;
*** TRANSPOSE ON AND OFF TREATMENT LSM                           ***;
********************************************************************;

proc transpose data = lsm_mmrm2
               out  = lsm_mmrm2_wide (drop=_name_)
               prefix = lsm;
  by rdrate sim_run visitn groupn;
  var estimate;
  id disc;
run;


********************************************************************;
*** CARCULATE PROPORTIONS FOR SUCCESSFULLY FITTED SIMS           ***;
********************************************************************;

proc sort data = ds;
  by rdrate sim_run visitn groupn disc subjid;
run;

data lsm_mmrm2_om (keep = rdrate sim_run visitn groupn ntotal ndisc: prop:);
  set ds (where = (mmrm2=1));
  retain ntotal ndisc0 ndisc1 0;
  by rdrate sim_run visitn groupn disc subjid;

  if first.groupn then do;   *** RESET COUNT FOR EACH GROUP WITHIN EACH VISIT ***;
    ntotal = 0; 
    ndisc0 = 0;
    ndisc1 = 0;
  end;

  ntotal = ntotal + 1;                      *** COUNT ALL SUBJECTS ***;
  select(disc);
    when (0) ndisc0 = ndisc0 + 1;
    when (1) ndisc1 = ndisc1 + 1;
    otherwise;
  end;

  if last.groupn then do;    *** OUTPUT FINAL COUNTS ***;
    prop0 = ndisc0 / ntotal;
    prop1 = ndisc1 / ntotal;    
    output; 
  end;

run;


********************************************************************;
*** PICK OUT VC MATRIX FOR SUCCESSFULLY FITTED SIMS              ***;
********************************************************************;

data lsm_mmrm2_vc;
  set lsm_mmrm2;
  by rdrate sim_run visitn groupn;
  retain row lsm_vc00 lsm_vc01 lsm_vc11 .;

  array covs[20] cov1-cov20;

  if first.sim_run then row = 0;
  if first.groupn then do;
    lsm_vc00 = .;
	lsm_vc01 = .;
	lsm_vc11 = .;
  end;

  row = row + 1;

  *** PICK OUT ON AND OFF TREATMENT VARIANCES AND COVARIANCE ***;
  if disc = 0 then lsm_vc00 = covs[row];
  else if disc = 1 then do;
    lsm_vc11 = covs[row];
    lsm_vc01 = covs[row-1];
  end;

  if last.groupn then output;
  keep rdrate sim_run visitn groupn lsm_vc00 lsm_vc01 lsm_vc11;
 
run;


********************************************************************;
*** MERGE PROPORTIONS LSMS AND VCS AND CALCULATE POLICY          ***;
********************************************************************;

data lsm_mmrm2_policy;
  merge lsm_mmrm2_om
        lsm_mmrm2_wide
        lsm_mmrm2_vc;
  by rdrate sim_run visitn groupn;

  lsm_policy_est = prop0*lsm0 + prop1*lsm1;
  
  lsm_policy_var   = (lsm_vc00*(prop0**2)) + (2*prop0*prop1*lsm_vc01) + (lsm_vc11*(prop1**2));
  lsm_policy_se    = sqrt( lsm_policy_var );
  lsm_policy_lower = lsm_policy_est - probit(0.975)*lsm_policy_se;
  lsm_policy_upper = lsm_policy_est + probit(0.975)*lsm_policy_se;

  lsm_policy_var_adj   = lsm_policy_var + ((1/ntotal)*prop0*prop1*((lsm1-lsm0)**2));
  lsm_policy_se_adj    = sqrt( lsm_policy_var_adj );
  lsm_policy_lower_adj = lsm_policy_est - probit(0.975)*lsm_policy_se_adj;
  lsm_policy_upper_adj = lsm_policy_est + probit(0.975)*lsm_policy_se_adj;

run;


********************************************************************;
*** CALC DIFS IN LSM                                             ***;
********************************************************************;

data dif_mmrm2_policy;
  merge lsm_mmrm2_policy (where=(groupn=1) rename=(ntotal = n1 lsm_policy_est=lsm_policy1 lsm_policy_var=lsm_policy_var1 lsm_policy_var_adj=lsm_policy_var_adj1))
        lsm_mmrm2_policy (where=(groupn=2) rename=(ntotal = n2 lsm_policy_est=lsm_policy2 lsm_policy_var=lsm_policy_var2 lsm_policy_var_adj=lsm_policy_var_adj2));
  by rdrate sim_run visitn;

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

  keep rdrate sim_run visitn _visitn groupn _groupn n1 n2 dif_policy_:;

run;


********************************************************************;
*** FIT MMRM1 TO REMAINING DATA WITH CORRECT BASELINE MEAN       ***;
********************************************************************;

%macro mmrm1(rdrate=);

  %ods_off();

  *** SUBSET TO MMRM1 DATA AND COUNT SIMS ***;
  data ds_mmrm1;
    set ds (where = (rdrate = &rdrate. and mmrm2=0)) end = eof;
    retain sim_count 0;
    by rdrate sim_run;
    if first.sim_run then sim_count + 1;
	if eof then call symput("nsims", put(sim_count,best.));
  run;

  %do ii = 1 %to &nsims.;

  *** SQL TO CREATE MACRO VARIABLE FOR BASELINE MEAN ***;
  proc sql noprint;
    select mean(baseline_var) into :mbl
    from ds_mmrm1
    where sim_count = &ii. and visitn = 1;
  quit;

  *** MMRM1 WITH LSM AT CORRECT BASELINE VALUE ***;
  proc mixed data = ds_mmrm1 (where = (sim_count = &ii.));
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

data lsm_mmrm1_policy;
  set lsm_r1
      lsm_r2
      lsm_r3
      lsm_r4
      lsm_r5
      lsm_r6;
  by rdrate sim_run visitn groupn;

  lsm_policy_est   = estimate;

  lsm_policy_var   = stderr*2;
  lsm_policy_se    = stderr;
  lsm_policy_lower = lower;
  lsm_policy_upper = upper;

  lsm_policy_var_adj   = stderr*2;
  lsm_policy_se_adj    = stderr;
  lsm_policy_lower_adj = lower;
  lsm_policy_upper_adj = upper;

  keep rdrate sim_run visitn groupn lsm_:;

run; 

data dif_mmrm1_policy;
  set dif_r1
      dif_r2
      dif_r3
      dif_r4
      dif_r5
      dif_r6;
  by rdrate sim_run visitn groupn;

  dif_policy_est   = estimate;

  dif_policy_var   = stderr*2;
  dif_policy_se    = stderr;
  dif_policy_lower = lower;
  dif_policy_upper = upper;

  dif_policy_var_adj   = stderr*2;
  dif_policy_se_adj    = stderr;
  dif_policy_lower_adj = lower;
  dif_policy_upper_adj = upper;

  keep rdrate sim_run visitn _visitn groupn _groupn dif_:;

run; 

%ods_off();
proc datasets lib = work nolist;
  delete ds_mmrm1 lsm_sc: lsm_r: dif_sc: dif_r:;
quit;
%ods_on();


********************************************************************;
*** STACK RESULTS FROM MMRM2 AND MMRM1 MODELS INTO FINAL RESULTS ***;
********************************************************************;

data lsm_policy;
  set lsm_mmrm2_policy (in = in1)
      lsm_mmrm1_policy (in = in2);
  by  rdrate sim_run visitn groupn;

  length fit $50;
  if in1 then do;
    fitn = 1;
	fit  = "Full MMRM2";
  end;
  else if in2 then do;
    fitn = 2;
	fit  = "Reduced to MMRM1";
  end;

run;


data dif_policy;
  set dif_mmrm2_policy (in = in1)
      dif_mmrm1_policy (in = in2);
  by  rdrate sim_run visitn groupn;

  length fit $50;
  if in1 then do;
    fitn = 1;
	fit  = "Full MMRM2";
  end;
  else if in2 then do;
    fitn = 2;
	fit  = "Reduced to MMRM1";
  end;

run;


********************************************************************;
*** ORGANISE RESULTS AND OUTPUT TO PERM LOCATION                 ***;
********************************************************************;

proc sort data = lsm_policy (keep = rdrate sim_run groupn visitn fitn fit lsm_policy_:)
          out  = results.&scenario._mmrm2_lsm_part&part.;
  by rdrate sim_run groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn fitn fit dif_policy_:)
          out  = results.&scenario._mmrm2_dif_part&part.;
  by rdrate sim_run visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY                                         ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: lsm: dif: prop: om: fit: ;
quit;














































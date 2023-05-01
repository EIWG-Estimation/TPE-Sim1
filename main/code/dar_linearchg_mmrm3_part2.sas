/*******************************************************************
| Name       : dar_linearchg_mmrm3_part2.sas
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

%let scenario = dar_linearchg;
%let part     = 2;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";


********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

proc sql noprint;
  create table ds as
    select *, 1-ontrt as disc 
    from data.&scenario.
    where 5000 lt sim_run le 10000
    order by rdrate, sim_run, groupn, subjid, visitn, disc;
quit;


********************************************************************;
*** TRANSPOSE TO WIDE FORMAT TO CODE PATTERNS                    ***;
********************************************************************;

data ds1;
  set ds;
  by rdrate sim_run groupn subjid;
  retain row 0 y1-y5 d1-d5 w1-w5 .;   *** RETAIN VARIABLES OVER ROWS ***;

  array y[5] y1-y5;   *** ARRAY TO HOLD RESPONSES ***;
  array d[5] d1-d5;   *** ARRAY TO HOLD DISCONTINUATION INDICATORS ***;
  array w[5] w1-w5;   *** ARRAY TO HOLD WITHDRAWAL INDICATORS ***;
  array p[5] p1-p5;   *** ARRAY TO HOLD PATTERN AT EACH VISIT ***;

  if first.subjid then do;   *** FOR EACH SUBJECT RESET THE ARRAYS ***;
    row = 0;
	do j = 1 to 5;
     y[j] = .;
	 d[j] = .;
	 w[j] = .;
	 p[j] = .;
	end;
  end;

  row    = row + 1;     *** COUNT THE ROW ***;
  y[row] = response;    
  d[row] = disc;        *** DISCONTINUED FLAG ***;
  w[row] = is_missing;  *** WITHDRAWAL FLAG ***;

  if last.subjid then do;   *** ONLY OUTPUT FINAL ROW PER PATIENT ***;
  
    disctime = 5 - sum(of d1-d5);
    withtime = 5 - sum(of w1-w5);
  
    do j = 1 to 5;
      p[j] = ifn(j <= disctime, 6, disctime+1); *** CREATE PATTERN AT EACH VISIT (LEVEL 6 = ON TREATMENT) ***;
    end; 
  
    output;
  
  end;
  keep rdrate sim_run subjid groupn group baseline_var disctime withtime y1-y5 d1-d5 w1-w5 p1-p5;

run;


********************************************************************;
*** COUNT THE DATA FOR EACH PATTERN AT EACH VISIT AND FLAG ZEROS ***;
********************************************************************;

data ds2;
  set ds1;
  by rdrate sim_run groupn;
  retain pv11-pv15 pv21-pv25 pv31-pv35 pv41-pv45 pv51-pv55 pv61-pv65 0;  

  array y[5] y1-y5;        *** SUBJECTS DATA AT EACH VISIT ***;
  array p[5] p1-p5;        *** SUBJECTS PATTERN AT EACH VISIT ***; 
  array pv[6,5] pv11-pv15
                pv21-pv25
                pv31-pv35
                pv41-pv45
                pv51-pv55
                pv61-pv65;  *** PATTERN BY VISIT DATA COUNTS ***;
    
  if first.groupn then do; 
    do i = 1 to 6;
    do j = 1 to 5;
      pv[i,j] = 0;   *** RESET ALL THE PATTERN VISIT FLAGS TO NO ***;
    end;
    end;
  end;    
  
  do j = 1 to 5;
    if y[j] ne . then pv[p[j], j] = pv[p[j], j] + 1;  *** COUNT OBSERVED DATA IN EACH PATTERN AT VISIT J ***;
  end;
  
  if last.groupn then do;   *** CHECK COUNTS OF EACH DISC PATTERN AT EACH VISIT IN EACH GROUPN FOR ISSUES ***;
  
    p_i1 = min(of pv11-pv15)=0;  *** ANY ZERO COUNTS IN PV11-PV15 THEN ISSUE WITH NO PAT 1 DATA AT A VISIT ***;
    p_i2 = min(of pv22-pv25)=0;  *** ANY ZERO COUNTS IN PV22-PV25 THEN ISSUE WITH NO PAT 2 DATA AT A VISIT ***;
    p_i3 = min(of pv33-pv35)=0;  *** ANY ZERO COUNTS IN PV33-PV35 THEN ISSUE WITH NO PAT 3 DATA AT A VISIT ***;
    p_i4 = min(of pv44-pv45)=0;  *** ANY ZERO COUNTS IN PV44-PV45 THEN ISSUE WITH NO PAT 4 DATA AT A VISIT ***;
    p_i5 = min(of pv55-pv55)=0;  *** ANY ZERO COUNTS IN PV55 ONLY THEN ISSUE WITH NO PAT 5 DATA AT A VISIT ***;
    
    d_i1 = sum(pv11) = 0;                     *** ZERO PV11 DATA THEN ISSUE WITH NO DISC DATA AT VISIT 1 ***;
    d_i2 = sum(pv12,pv22) = 0;                *** ZERO PV12-PV22 DATA THEN ISSUE WITH NO DISC DATA AT VISIT 2 ***;
    d_i3 = sum(pv13,pv23,pv33) = 0;           *** ZERO PV13-PV33 DATA THEN ISSUE WITH NO DISC DATA AT VISIT 3 ***;
    d_i4 = sum(pv14,pv24,pv34,pv44) = 0;      *** ZERO PV14-PV44 DATA THEN ISSUE WITH NO DISC DATA AT VISIT 4 ***;
    d_i5 = sum(pv15,pv25,pv35,pv45,pv55) = 0; *** ZERO PV15-PV55 DATA THEN ISSUE WITH NO DISC DATA AT VISIT 5 ***;
        
    output;
  end;
  
  keep rdrate sim_run groupn p_i: d_i:;

run;


********************************************************************;
*** CHECK FOR ISSUE IN EACH GROUP FOR EACH PATTERN AT EACH VISIT ***;
********************************************************************;

proc sql noprint;

   create table ds3 as
     select rdrate, sim_run, 
            max(p_i1) as p_i1, max(p_i2) as p_i2, max(p_i3) as p_i3, max(p_i4) as p_i4, max(p_i5) as p_i5,    
            max(d_i1) as d_i1, max(d_i2) as d_i2, max(d_i3) as d_i3, max(d_i4) as d_i4, max(d_i5) as d_i5    
     from ds2
     group by rdrate, sim_run;

quit;


********************************************************************;
*** CHANGE BACK TO LONG FORMAT FOR FITTING MMRM MODELS           ***;
********************************************************************;

data ds4;
  merge ds1
        ds3;
  by rdrate sim_run;
  
  array y[5] y1-y5;   *** ARRAY TO HOLD RESPONSES ***;
  array d[5] d1-d5;   *** ARRAY TO HOLD DISCONTINUATION INDICATORS ***;
  array p[5] p1-p5;   *** ARRAY TO HOLD PATTERNS AT EACH VISIT ***;
  array w[5] w1-w5;   *** ARRAY TO HOLD WITHDRAWAL INDICATORS ***;

  do j = 1 to 5;
    visitn   = j;
    response = y[j];
    disc     = d[j];
    pattern  = p[j];
    with     = w[j];
    if response ne . then change = response - baseline_var;
    else change = .;
    output;
  end;

  keep rdrate sim_run subjid groupn group visitn baseline_var 
       response change disc disctime pattern with withtime p_i1-p_i5 d_i1-d_i5;

run;


********************************************************************;
*** COMBINE PATTERNS WITH ISSUES                                 ***;
********************************************************************;

data ds5;
  set ds4;
  
  issues = cat(p_i1, p_i2, p_i3, p_i4, p_i5); *** 32 POSSIBLE PATTERN ISSUES ***;     
  select (issues);
  
    when ("11111") do; if pattern in (1 2 3 4 5 6) then patcode = 123456; end; *** ONLY 1 PATTERN POSSIBLE (MMRM1) ***;
  
    when ("11110") do; if pattern in (1 2 3 4 5) then patcode = 12345; else patcode = pattern; end;  *** ONLY 2 PATTERNS POSSIBLE (MMRM2) ***;
    when ("11101") do; if pattern in (1 2 3 4 5) then patcode = 12345; else patcode = pattern; end; 
    when ("11011") do; if pattern in (1 2 3 4 5) then patcode = 12345; else patcode = pattern; end; 
    when ("10111") do; if pattern in (1 2 3 4 5) then patcode = 12345; else patcode = pattern; end;
    when ("01111") do; if pattern in (1 2 3 4 5) then patcode = 12345; else patcode = pattern; end;
    
    when ("11010") do; if pattern in (1 2 3 4) then patcode = 1234; else patcode = pattern; end; *** 3 PATTERNS POSSIBLE ***;
    when ("10110") do; if pattern in (1 2 3 4) then patcode = 1234; else patcode = pattern; end;
    when ("11100") do; if pattern in (1 2 3 4) then patcode = 1234; else patcode = pattern; end; 
    when ("01110") do; if pattern in (1 2 3 4) then patcode = 1234; else patcode = pattern; end;
    when ("00111") do; if pattern in (2 3 4 5) then patcode = 2345; else patcode = pattern; end;

    when ("11001") do; if pattern in (1 2 3) then patcode = 123; else if pattern in (4 5) then patcode = 45; else patcode = pattern; end;
    when ("10101") do; if pattern in (1 2 3) then patcode = 123; else if pattern in (4 5) then patcode = 45; else patcode = pattern; end;
    when ("01101") do; if pattern in (1 2 3) then patcode = 123; else if pattern in (4 5) then patcode = 45; else patcode = pattern; end;
    when ("10011") do; if pattern in (1 2) then patcode = 12; else if pattern in (3 4 5) then patcode = 345; else patcode = pattern; end;
    when ("01011") do; if pattern in (1 2) then patcode = 12; else if pattern in (3 4 5) then patcode = 345; else patcode = pattern; end;
    
    when ("10100") do; if pattern in (1 2 3) then patcode = 123; else patcode = pattern; end;  *** 4 PATTERNS POSSIBLE ***;
    when ("11000") do; if pattern in (1 2 3) then patcode = 123; else patcode = pattern; end;
    when ("01100") do; if pattern in (1 2 3) then patcode = 123; else patcode = pattern; end;
    when ("00110") do; if pattern in (2 3 4) then patcode = 234; else patcode = pattern; end;
    when ("00011") do; if pattern in (3 4 5) then patcode = 345; else patcode = pattern; end;

    when ("10010") do; if pattern in (1 2) then patcode = 12; else if pattern in (3 4) then patcode = 34; else patcode = pattern; end;
    when ("01010") do; if pattern in (1 2) then patcode = 12; else if pattern in (3 4) then patcode = 34; else patcode = pattern; end;
    when ("10001") do; if pattern in (1 2) then patcode = 12; else if pattern in (4 5) then patcode = 45; else patcode = pattern; end;
    when ("01001") do; if pattern in (1 2) then patcode = 12; else if pattern in (4 5) then patcode = 45; else patcode = pattern; end;
    when ("00101") do; if pattern in (2 3) then patcode = 23; else if pattern in (4 5) then patcode = 45; else patcode = pattern; end;

    when ("10000") do; if pattern in (1 2) then patcode = 12; else patcode = pattern; end;  *** 5 PATTERNS POSSIBLE ***;
    when ("01000") do; if pattern in (1 2) then patcode = 12; else patcode = pattern; end;
    when ("00100") do; if pattern in (2 3) then patcode = 23; else patcode = pattern; end;
    when ("00010") do; if pattern in (3 4) then patcode = 34; else patcode = pattern; end;
    when ("00001") do; if pattern in (4 5) then patcode = 45; else patcode = pattern; end;
  
    when ("00000") do; patcode = pattern; end;  *** FULL 6 PATTERNS POSSIBLE (MMRM3) ***;

    otherwise;
 
  end;

  select(issues);
    when ("00000") fitn = 1;
    when ("10000", "01000", "00100", "00010", "00001") fitn = 2;
    when ("10100", "11000", "01100", "00110", "00011", "10010", "01010", "10001", "01001", "00101") fitn = 3;    
    when ("11010", "10110", "11100", "01110", "00111", "11001", "10101", "01101", "10011", "01011") fitn = 4;
    when ("11110", "11101", "11011", "10111", "01111") fitn = 5; 
    when ("11111") fitn = 6;
    otherwise;
  end;
    
  length fit $30;
  select(fitn);
    when (1) fit = "6 Patterns - MMRM3";
    when (2) fit = "5 Patterns";
    when (3) fit = "4 Patterns";
    when (4) fit = "3 Patterns";
    when (5) fit = "2 Patterns - MMRM2";
    when (6) fit = "1 Pattern - MMRM1";
    otherwise;
  end;
  
  drop p_i: d_i:;
  
run;


proc sort data = ds5
          out  = ds;
    by rdrate sim_run visitn groupn pattern subjid;
run;


********************************************************************;
*** CLEAN UP WORK LIBRARY                                        ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds1 ds2 ds3 ds4 ds5;
quit;


********************************************************************;
*** FIT MMRM3 USING PATTERN COVARIATE AND CORRECT BASELINE MEAN  ***;
********************************************************************;

%macro mmrm3(rdrate=);

  %ods_off(notes=N);
   
  %let nsims = 0;    *** IN CASE THERE ARE NO CASES IN DATA ***;
   
  *** COUNT SIMS ***;
  data ds_mmrm3;
    set ds (where = (rdrate = &rdrate.)) end = eof;
    retain sim_count 0;
    by rdrate sim_run fitn fit;
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
    by rdrate sim_run fit fitn;
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
    by rdrate sim_run fitn fit visitn groupn patcode;
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
  by rdrate sim_run fitn fit visitn groupn patcode;
run; 

%ods_off();
proc datasets lib = work nolist;
  delete lsm_sc: lsm_r: ds_mmrm3;
quit;
%ods_on();


********************************************************************;
*** CALCULATE PROPORTIONS BY PATTERN AND TRANSPOSE               ***;
********************************************************************;

proc freq data = ds noprint;
  by rdrate sim_run fitn fit visitn groupn;
  tables patcode / out = mmrm3_om;
run;

data lsm_mmrm3_om;
  set mmrm3_om;
  by rdrate sim_run fitn fit visitn groupn patcode;
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

  keep rdrate sim_run fitn fit visitn groupn n_pat:;
 
run;


********************************************************************;
*** TRANSPOSE THE LSM AND VC ESTIMATES                           ***;
********************************************************************;

data lsm_mmrm3_wide;
  set lsm_mmrm3;
  by rdrate sim_run fitn fit visitn groupn patcode;
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
  keep rdrate sim_run fitn fit visitn groupn lsm_pat: lsm_vc:;
 
run;


********************************************************************;
*** MERGE PROPORTIONS LSMS AND VCS AND CALCULATE POLICY          ***;
********************************************************************;

data lsm_policy;
  merge lsm_mmrm3_om
        lsm_mmrm3_wide;
  by rdrate sim_run fitn fit visitn groupn;
  
  ntotal = sum(of n_pat1-n_pat6);
  npats  = 6-nmiss(lsm_pat1,lsm_pat2,lsm_pat3,lsm_pat4,lsm_pat5,lsm_pat6);

  if npats = 1 then do;         *** FORMULA FOR POLICY AND VARIANCE ADJUSTMENT BASED ON PATTERNS ***;
  
    lsm_policy_est     = lsm_pat1;  *** IF ONLY ONE PAT POSSIBLE NO WEIGHING NEEDED ***;
    lsm_policy_var     = lsm_vc11;  
    lsm_policy_var_adj = lsm_vc11;  *** IF ONLY ONE PAT POSSIBLE NO DELTA METHOD ADJUSTMENT NEEDED ***;
    
  end;
  else if npats = 2 then do;
    
    p_pat1 = n_pat1 / ntotal;
    p_pat2 = n_pat2 / ntotal;  *** PAT2 IS ON TREATMENT ***;

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
    p_pat3 = n_pat3 / ntotal;  *** PAT3 IS ON TREATMENT ***;

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
    p_pat4 = n_pat4 / ntotal;  *** PAT4 IS ON TREATMENT ***;

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
    p_pat5 = n_pat5 / ntotal;  *** PAT5 IS ON TREATMENT ***;

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
    p_pat6 = n_pat6 / ntotal;  *** PAT6 IS ON TREATMENT ***;

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

  keep rdrate sim_run fitn fit visitn groupn ntotal lsm_policy_:;
  
run;


********************************************************************;
*** CALC DIFS IN LSM                                             ***;
********************************************************************;

data dif_policy;
  merge lsm_policy (where=(groupn=1) rename=(ntotal = n1 lsm_policy_est=lsm_policy1 lsm_policy_var=lsm_policy_var1 lsm_policy_var_adj=lsm_policy_var_adj1))
        lsm_policy (where=(groupn=2) rename=(ntotal = n2 lsm_policy_est=lsm_policy2 lsm_policy_var=lsm_policy_var2 lsm_policy_var_adj=lsm_policy_var_adj2));
  by rdrate sim_run fitn fit visitn;

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

  keep rdrate sim_run fitn fit visitn _visitn groupn _groupn n1 n2 dif_policy_:;

run;


********************************************************************;
*** ORGANISE RESULTS AND OUTPUT TO PERM LOCATION                 ***;
********************************************************************;

proc sort data = lsm_policy (keep = rdrate sim_run groupn visitn fitn fit lsm_policy_:)
          out  = results.&scenario._mmrm3_lsm_part&part.;
  by rdrate sim_run fitn fit groupn visitn;
run;

proc sort data = dif_policy (keep = rdrate sim_run groupn _groupn visitn _visitn fitn fit dif_policy_:)
          out  = results.&scenario._mmrm3_dif_part&part.;
  by rdrate sim_run fitn fit visitn;
run;


********************************************************************;
*** TIDY UP WORK LIBRARY                                         ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds: lsm: dif: om: mmrm: ;
quit;

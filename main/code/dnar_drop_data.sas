/*******************************************************************
| Name       : dnar_drop_data.sas
| Purpose    : Code to create the disc and pattern status 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 28OCT23 
|-------------------------------------------------------------------
| Notes:
|
|
*******************************************************************/;

********************************************************************;
*** SET UP ENVIRONMENT                                           ***;
********************************************************************;

%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";
%let scenario = dnar_drop;


********************************************************************;
*** READ IN DATA                                                 ***;
********************************************************************;

proc sql noprint;
  create table ds as
    select *, 1-ontrt as disc 
    from data.&scenario.
    where 0 lt sim_run le 10000 
    order by rdrate, sim_run, groupn, subjid, visitn, disc;
quit;


********************************************************************;
*** TRANSPOSE TO WIDE FORMAT TO CODE PATTERNS                    ***;
********************************************************************;

data ds1;
  set ds;
  by rdrate sim_run groupn subjid;
  retain row 0 yf1-yf5 y1-y5 d1-d5 w1-w5 .;   *** RETAIN VARIABLES OVER ROWS ***;

  array yf[5] yf1-yf5;  *** ARRAY TO HOLD FULL RESPONSES ***;

  array y[5] y1-y5;     *** ARRAY TO HOLD OBSERVED RESPONSES ***;
  array d[5] d1-d5;     *** ARRAY TO HOLD DISCONTINUATION INDICATORS ***;
  array w[5] w1-w5;     *** ARRAY TO HOLD WITHDRAWAL INDICATORS ***;
  array p[5] p1-p5;     *** ARRAY TO HOLD PATTERN AT EACH VISIT ***;

  if first.subjid then do;   *** FOR EACH SUBJECT RESET THE ARRAYS ***;
    row = 0;
	do j = 1 to 5;
	 yf[j] = .;
     y[j] = .;
	 d[j] = .;
	 w[j] = .;
	 p[j] = .;
	end;
  end;

  row    = row + 1;     *** COUNT THE ROW ***;
  yf[row] = response_full;
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
  keep rdrate sim_run subjid groupn group baseline_var disctime withtime yf1-yf5 y1-y5 d1-d5 w1-w5 p1-p5;

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

    array d_i[5] d_i1 d_i2 d_i3 d_i4 d_i5;               *** NO DISC DATA FLAGS ***;
    array e_i[5] e_i1 e_i2 e_i3 e_i4 e_i5 (0 0 0 0 0);   *** ESTIMATION ISSUE FLAGS ***;

    do i = 1 to 5;
      e_i[i] = 0;
      do j = i to 5;
        if pv[i,j] = 0 and d_i[j] = 1 then e_i[i] = 1;   *** NO PATTERN DATA AT VISIT AND NO OTHER OFF TREATMENT DATA AT VISIT TO USE ***; 
      end;
    end;
        
    output;
  end;
  
run;


********************************************************************;
*** CHECK FOR ISSUE IN EACH GROUP FOR EACH PATTERN AT EACH VISIT ***;
********************************************************************;

proc sql noprint;

   create table ds3 as
     select rdrate, sim_run, 
            max(p_i1) as p_i1, max(p_i2) as p_i2, max(p_i3) as p_i3, max(p_i4) as p_i4, max(p_i5) as p_i5,    
            max(d_i1) as d_i1, max(d_i2) as d_i2, max(d_i3) as d_i3, max(d_i4) as d_i4, max(d_i5) as d_i5,
            max(e_i1) as e_i1, max(e_i2) as e_i2, max(e_i3) as e_i3, max(e_i4) as e_i4, max(e_i5) as e_i5            
     from ds2
     group by rdrate, sim_run;

quit;


********************************************************************;
*** COMBINE PATTERNS WITH ISSUES                                 ***;
********************************************************************;

data ds;
  merge ds1
        ds3;
  by rdrate sim_run;
   
  issues = cat(p_i1, p_i2, p_i3, p_i4, p_i5); *** 32 POSSIBLE PATTERN ISSUES ***;     
     
  array pattern[5] p1-p5; 
  array disccode[5] disc1-disc5;     
  array patcode[5] pat1-pat5;
       
  *** CODE PATTERNS INTO DISC CODES FOR MMRM2 AND MI2 MODELS BASED ON IDENTIFIABILITY ***;     
  select (issues);      
    when ("11111") do i = 1 to 5; if pattern[i] in (1 2 3 4 5 6) then disccode[i] = 123456; end; *** ONLY 1 PATTERN POSSIBLE ***;
    otherwise do i = 1 to 5; if pattern[i] in (1 2 3 4 5) then disccode[i] = 12345; else disccode[i] = pattern[i]; end; *** 2 PATTERNS POSSIBLE ***;  
  end;
 
  *** CODE PATTERNS INTO PAT CODES FOR MMRM3 AND MI3 MODELS BASED ON IDENTIFIABILITY ***;
  select (issues); 
  
    when ("11111") do i = 1 to 5; if pattern[i] in (1 2 3 4 5 6) then patcode[i] = 123456; end; *** ONLY 1 PATTERN POSSIBLE ***;
  
    when ("11110") do i = 1 to 5; if pattern[i] in (1 2 3 4 5) then patcode[i] = 12345; else patcode[i] = pattern[i]; end;  *** ONLY 2 PATTERNS POSSIBLE ***;
    when ("11101") do i = 1 to 5; if pattern[i] in (1 2 3 4 5) then patcode[i] = 12345; else patcode[i] = pattern[i]; end; 
    when ("11011") do i = 1 to 5; if pattern[i] in (1 2 3 4 5) then patcode[i] = 12345; else patcode[i] = pattern[i]; end; 
    when ("10111") do i = 1 to 5; if pattern[i] in (1 2 3 4 5) then patcode[i] = 12345; else patcode[i] = pattern[i]; end;
    when ("01111") do i = 1 to 5; if pattern[i] in (1 2 3 4 5) then patcode[i] = 12345; else patcode[i] = pattern[i]; end;
    
    when ("11010") do i = 1 to 5; if pattern[i] in (1 2 3 4) then patcode[i] = 1234; else patcode[i] = pattern[i]; end; *** 3 PATTERNS POSSIBLE ***;
    when ("10110") do i = 1 to 5; if pattern[i] in (1 2 3 4) then patcode[i] = 1234; else patcode[i] = pattern[i]; end;
    when ("11100") do i = 1 to 5; if pattern[i] in (1 2 3 4) then patcode[i] = 1234; else patcode[i] = pattern[i]; end; 
    when ("01110") do i = 1 to 5; if pattern[i] in (1 2 3 4) then patcode[i] = 1234; else patcode[i] = pattern[i]; end;
    when ("00111") do i = 1 to 5; if pattern[i] in (2 3 4 5) then patcode[i] = 2345; else patcode[i] = pattern[i]; end;

    when ("11001") do i = 1 to 5; if pattern[i] in (1 2 3) then patcode[i]= 123; else if pattern[i] in (4 5)   then patcode[i] = 45;  else patcode[i] = pattern[i]; end;
    when ("10101") do i = 1 to 5; if pattern[i] in (1 2 3) then patcode[i]= 123; else if pattern[i] in (4 5)   then patcode[i] = 45;  else patcode[i] = pattern[i]; end;
    when ("01101") do i = 1 to 5; if pattern[i] in (1 2 3) then patcode[i]= 123; else if pattern[i] in (4 5)   then patcode[i] = 45;  else patcode[i] = pattern[i]; end;
    when ("10011") do i = 1 to 5; if pattern[i] in (1 2)   then patcode[i]= 12;  else if pattern[i] in (3 4 5) then patcode[i] = 345; else patcode[i] = pattern[i]; end;
    when ("01011") do i = 1 to 5; if pattern[i] in (1 2)   then patcode[i]= 12;  else if pattern[i] in (3 4 5) then patcode[i] = 345; else patcode[i] = pattern[i]; end;
    
    when ("10100") do i = 1 to 5; if pattern[i] in (1 2 3) then patcode[i] = 123; else patcode[i] = pattern[i]; end;  *** 4 PATTERNS POSSIBLE ***;
    when ("11000") do i = 1 to 5; if pattern[i] in (1 2 3) then patcode[i] = 123; else patcode[i] = pattern[i]; end;
    when ("01100") do i = 1 to 5; if pattern[i] in (1 2 3) then patcode[i] = 123; else patcode[i] = pattern[i]; end;
    when ("00110") do i = 1 to 5; if pattern[i] in (2 3 4) then patcode[i] = 234; else patcode[i] = pattern[i]; end;
    when ("00011") do i = 1 to 5; if pattern[i] in (3 4 5) then patcode[i] = 345; else patcode[i] = pattern[i]; end;

    when ("10010") do i = 1 to 5; if pattern[i] in (1 2) then patcode[i] = 12; else if pattern[i] in (3 4) then patcode[i] = 34; else patcode[i] = pattern[i]; end;
    when ("01010") do i = 1 to 5; if pattern[i] in (1 2) then patcode[i] = 12; else if pattern[i] in (3 4) then patcode[i] = 34; else patcode[i] = pattern[i]; end;
    when ("10001") do i = 1 to 5; if pattern[i] in (1 2) then patcode[i] = 12; else if pattern[i] in (4 5) then patcode[i] = 45; else patcode[i] = pattern[i]; end;
    when ("01001") do i = 1 to 5; if pattern[i] in (1 2) then patcode[i] = 12; else if pattern[i] in (4 5) then patcode[i] = 45; else patcode[i] = pattern[i]; end;
    when ("00101") do i = 1 to 5; if pattern[i] in (2 3) then patcode[i] = 23; else if pattern[i] in (4 5) then patcode[i] = 45; else patcode[i] = pattern[i]; end;

    when ("10000") do i = 1 to 5; if pattern[i] in (1 2) then patcode[i] = 12; else patcode[i] = pattern[i]; end;  *** 5 PATTERNS POSSIBLE ***;
    when ("01000") do i = 1 to 5; if pattern[i] in (1 2) then patcode[i] = 12; else patcode[i] = pattern[i]; end;
    when ("00100") do i = 1 to 5; if pattern[i] in (2 3) then patcode[i] = 23; else patcode[i] = pattern[i]; end;
    when ("00010") do i = 1 to 5; if pattern[i] in (3 4) then patcode[i] = 34; else patcode[i] = pattern[i]; end;
    when ("00001") do i = 1 to 5; if pattern[i] in (4 5) then patcode[i] = 45; else patcode[i] = pattern[i]; end;
  
    when ("00000") do i = 1 to 5; patcode[i] = pattern[i]; end;  *** FULL 6 PATTERNS POSSIBLE (MMRM3) ***;

    otherwise;
 
  end;

  *** FIX ESTIMATION ISSUES FOR 2 AND 3 MODELS THAT HAVE PATTERN ISSUE BUT NO AVAILABLE DISC PATTERN DATA YET E.G VISIT 1 ***;
  if issues ne "11111" then do;              
    if p_i1 = 1 and e_i1 = 1 then do; disc1 = 6; pat1 = 6; end;
    if p_i2 = 1 and e_i2 = 1 then do; disc2 = 6; pat2 = 6; end;
    if p_i3 = 1 and e_i3 = 1 then do; disc3 = 6; pat3 = 6; end;
    if p_i4 = 1 and e_i4 = 1 then do; disc4 = 6; pat4 = 6; end;
    if p_i5 = 1 and e_i5 = 1 then do; disc5 = 6; pat5 = 6; end;
  end;  
  
  *** CREATE VARIABLE SHOWING HOW MANY PATTERNS POSSIBLE ***;
  select(issues);
    when ("11111") patmaxn = 1;
    when ("11110", "11101", "11011", "10111", "01111") patmaxn = 2; 
    when ("11010", "10110", "11100", "01110", "00111", "11001", "10101", "01101", "10011", "01011") patmaxn = 3;
    when ("10100", "11000", "01100", "00110", "00011", "10010", "01010", "10001", "01001", "00101") patmaxn = 4;    
    when ("10000", "01000", "00100", "00010", "00001") patmaxn = 5;
    when ("00000") patmaxn = 6;
    otherwise;
  end;
    
  length patmax $30;
  select(patmaxn);
    when (1) patmax = "1 Pattern";
    when (2) patmax = "2 Patterns";
    when (3) patmax = "3 Patterns";
    when (4) patmax = "4 Patterns";
    when (5) patmax = "5 Patterns";
    when (6) patmax = "6 Patterns";
    otherwise;
  end;
 
run;   

*** CHECKS ***;
/* proc freq data = ds; */
/*   tables (disc1 pat1)*rdrate / norow nocol; */
/*   tables (disc2 pat2)*rdrate / norow nocol; */
/*   tables (disc3 pat3)*rdrate / norow nocol; */
/*   tables (disc4 pat4)*rdrate / norow nocol; */
/*   tables (disc5 pat5)*rdrate / norow nocol; */
/* run; */
/*      */
/* proc sort data=ds out=ds4 nodupkey; */
/*   by rdrate sim_run; */
/* run; */
/*  */
/* proc freq data = ds4; */
/*   tables patmax*rdrate / norow nocol; */
/* run; */
 
********************************************************************;
*** CREATE PERMANENT DATASET                                     ***;
********************************************************************;

data data.&scenario._data;
  set ds;
  keep rdrate sim_run subjid groupn group baseline_var	
       yf1-yf5 y1-y5 
       disctime d1-d5 p1-p5 
       withtime w1-w5
       issues p_i1-p_i5	d_i1-d_i5 e_i1-e_i5 
       patmaxn patmax	
       disc1-disc5	pat1-pat5;	
run;


********************************************************************;
*** CLEAN UP WORK LIBRARY                                        ***;
********************************************************************;

proc datasets lib = work nolist;
  delete ds:;
quit;




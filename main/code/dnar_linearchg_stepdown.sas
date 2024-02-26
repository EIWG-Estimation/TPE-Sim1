/*******************************************************************
| Name       : dnar_linearchg_stepdown.sas
| Purpose    : Code to calculate number of step downs 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 05AUG23 
|-------------------------------------------------------------------
| Notes:
| 
*******************************************************************/;

********************************************************************;
*** SECTION 0: SET UP ENVIRONMENT                                ***;
********************************************************************;

%let scenario = dnar_linearchg;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";


********************************************************************;
*** SECTION 1: READ AND STACK DATA                               ***;
********************************************************************;

data dif_all;
  set results.&scenario._mmrm2_dif_part1 (in = in2 where=(visitn = 5))
      results.&scenario._mmrm2_dif_part2 (in = in2 where=(visitn = 5))
      results.&scenario._mmrm3_dif_part1 (in = in3 where=(visitn = 5))
      results.&scenario._mmrm3_dif_part2 (in = in3 where=(visitn = 5));
  
  if in2 then approachn = 1;
  else if in3 then approachn = 2;
  select(approachn);
    when (1) approach= "MMRM2";
    when (2) approach= "MMRM3";
    otherwise;    
  end;
  
run;  
  
  
********************************************************************;
*** SECTION 2: GET FREQUENCIES AND PERCENTAGES                   ***;
********************************************************************;
  
proc freq data = dif_all;
  by approach;
  tables fit*rdrate;
run;
    
 
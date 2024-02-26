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
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/null/code/analysis_setup.sas";


********************************************************************;
*** SECTION 1: READ IN DATA                                      ***;
********************************************************************;

proc sort data = results.&scenario._summary_dif (where=(visitn = 5))
          out  = dif_summary;
  by methodn method rdrate;
run;  

data dif_summary;  
  set dif_summary;
  
  if method in ("MMRM2", "MMRM3") then dif_policy_type1_pct = 100*dif_policy_type1_adj_mean;
  else dif_policy_type1_pct = 100*dif_policy_type1_mean;

run;
  
********************************************************************;
*** SECTION 2: TRANSPOSE                                         ***;
********************************************************************;
  
proc transpose data   = dif_summary
               out    = dif_type1
               prefix = rdrate;  
  by methodn method dif_policy_type1_n;
  id rdrate;
  var dif_policy_type1_pct;
run;
  
  
  
proc print data = dif_type1 noobs label;
  var method dif_policy_type1_n rdrate1-rdrate6;
  label method = "Model"
        dif_policy_type1_n = "N";
run;
    
 
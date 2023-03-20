/*******************************************************************
| Name       : dar_drop_data.sas
| Purpose    : Code to read in .CSV simulation data into SAS format 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 10SEP21 
|-------------------------------------------------------------------
| Notes:
|
|
*******************************************************************/;

********************************************************************;
*** SECTION 0: SET UP ENVIRONMENT                                ***;
********************************************************************;

%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/code/analysis_setup.sas";
%let scenario = dar_drop;

*** SET UP FILE NAME TO THE CSV AND LIBRARY FOR STORING SAS DATA ***;
filename csv0 "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/data/&scenario._true.csv";

filename csv1 "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/data/&scenario._p01.csv";
filename csv2 "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/data/&scenario._p02.csv";
filename csv3 "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/data/&scenario._p03.csv";
filename csv4 "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/data/&scenario._p04.csv";
filename csv5 "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/data/&scenario._p05.csv";
filename csv6 "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/data/&scenario._p06.csv";

********************************************************************;
*** SECTION 1: MACRO TO READ IN DATA                             ***;
********************************************************************;

%macro read_csv(csvfile=, dsetname=);

*** READ IN DATA ***;
data work.dset1 (drop = _:);
  infile &csvfile. dsd missover dlm="," firstobs=2;
  input sim_run  
        subjid  
        visitn 
        mean_trt
        mean_ctl
        response_ontrt
        group $
        baseline_var
        ontrt
        response_offtrt
        is_missing
        _response $;

  *** HANDLE MISSING NAs ***;
  if strip(_response) = "NA" then response = .;
  else response = input(_response, best12.);

run;

*** SORT AND REMOVE VISIT ZERO ***;
proc sort data = work.dset1
          out  = work.dset2;
  by sim_run group subjid visitn;
  where visitn in (1 2 3 4 5);
run;

*** CREATE EXTRA OUTCOME VARIABLES ***;
data work.dset3;
  set work.dset2; 
  by sim_run group subjid visitn;

  *** CREATE UNIQUE SUBJECT ID ACROSS ALL SIMULATIONS ***;
  usubjid = 10000*sim_run + subjid; 

  *** CREATE NUMERIC VARIABLE FOR GROUP TO GET DIFFS TRT vs CTL ***;
  groupn = 2*(group="ctl") + 1*(group="trt");

  *** CREATE CHANGE FROM BASELINE VARIABLE IF RESPONSE AND BASELINE NON MISSING ***;
  if nmiss(response,baseline_var) = 0 then change = response - baseline_var;

  *** CREATE A VARIABLE WITH TRUE POLICY RESPONSE WITHOUT MISSINGNESS ***;
  response_full = (ontrt=1)*(response_ontrt) + (ontrt=0)*(response_offtrt);

  *** CREATE FULL CHANGE FROM BASELINE VARIABLE ***;
  change_full = response_full - baseline_var;

run;

*** CREATE DISCONTINUATION TIME VARIABLES ***;
data work.dset4 (keep = sim_run group subjid disctime);
  set work.dset3;
  by sim_run group subjid visitn;

  retain disctime 0;
  if first.subjid then disctime = 0;
  if ontrt = 1 then disctime = visitn;
  if last.subjid then output;

run;

*** REMERGE TO ORIGINAL DATA ***;
data work.dset5;
  merge dset3
        dset4;
  by sim_run group subjid;
run;  

*** OUTPUT WITH SPECIFIED NAME ***;
data &dsetname.;
  set work.dset5;
run;

*** CLEAN UP WORK ***;
proc datasets lib = work nolist;
  delete dset:;
quit;

%mend;


********************************************************************;
*** SECTION 2: CALL MACRO TO READ IN ALL CSVS                    ***;
********************************************************************;
%read_csv(csvfile=csv0, dsetname=data.&scenario._true);

%read_csv(csvfile=csv1, dsetname=&scenario._p01);
%read_csv(csvfile=csv2, dsetname=&scenario._p02);
%read_csv(csvfile=csv3, dsetname=&scenario._p03);
%read_csv(csvfile=csv4, dsetname=&scenario._p04);
%read_csv(csvfile=csv5, dsetname=&scenario._p05);
%read_csv(csvfile=csv6, dsetname=&scenario._p06);


********************************************************************;
*** SECTION 3: STACK INTO A SINGLE DATASET                       ***;
********************************************************************;

data data.&scenario.;
  set &scenario._p01 (in = in1)
      &scenario._p02 (in = in2)
      &scenario._p03 (in = in3)
      &scenario._p04 (in = in4)
      &scenario._p05 (in = in5)
      &scenario._p06 (in = in6);

  if in1 then rdrate = 1;
  else if in2 then rdrate = 2;
  else if in3 then rdrate = 3;
  else if in4 then rdrate = 4;
  else if in5 then rdrate = 5;
  else if in6 then rdrate = 6;

run;


********************************************************************;
*** SECTION 4: CLEAR UP WORK DATASET                             ***;
********************************************************************;

proc datasets lib = work nolist;
  delete &scenario._:;
quit;

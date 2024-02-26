/*******************************************************************
| Name       : analysis_setup.sas
| Purpose    : Code to set up libnames and call utility functions
|              and macros. 
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 15SEP21 
|-------------------------------------------------------------------
| Notes:
|
| 
*******************************************************************/;


*** SIMSTUDY LIBNAMES ***;
libname data    "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/null/data";
libname results "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/null/results"; 


*** MACRO TO TURN OFF ALL AUTOMATIC OUTPUT ***;
%macro ods_off(notes=N);
ods results off;
ods graphics off;
ods select none;
%if &notes.=N %then %do; options nonotes; %end;
%mend;


*** MACRO TO TURN ON ALL AUTOMATIC OUTPUT ***;
%macro ods_on();
ods results on;
ods graphics on;
ods select all;
options notes; 
%mend;

/*******************************************************************
| Name       : dar_drop_status.sas
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

%let scenario = dar_drop;
%inc "/hpawrk/tad66240/collaboration/eig_estimation/tpe-sim1/main/code/analysis_setup.sas";


********************************************************************;
*** SECTION 1: READ DATA                                         ***;
********************************************************************;

data ds1;
  set data.&scenario._data;
  where sim_run le 5000;
  
  array d[5] d1-d5;
  array w[5] w1-w5;
  
  array o[5];
  array x[5];
  array m[5];
  
  do j = 1 to 5;
    o[j] = d[j]=0 & w[j]=0;
    x[j] = d[j]=1 & w[j]=0;
    m[j] = w[j];
  end;
  keep groupn group rdrate sim_run o1-o5 x1-x5 m1-m5 d1-d5 w1-w5;
run;

proc means data = ds1 nway noprint;
  class rdrate groupn group;
  vars o1 x1 m1 
       o2 x2 m2 
       o3 x3 m3 
       o4 x4 m4  
       o5 x5 m5;
  output out = ss1 (where = (_stat_ = "MEAN"));
run;

data ss2 (keep = rdrate group: type: v:);
  set ss1 (in = in1 keep = groupn group rdrate o1-o5 rename = (o1=v1 o2=v2 o3=v3 o4=v4 o5=v5))
      ss1 (in = in2 keep = groupn group rdrate x1-x5 rename = (x1=v1 x2=v2 x3=v3 x4=v4 x5=v5))
      ss1 (in = in3 keep = groupn group rdrate m1-m5 rename = (m1=v1 m2=v2 m3=v3 m4=v4 m5=v5));
  by rdrate groupn group;
  
  if in1 then typen = 1;
  else if in2 then typen = 2;
  else if in3 then typen = 3;
  
  select (typen);
    when (1) type = "O";
    when (2) type = "X";
    when (3) type = "M";
    otherwise;
  end;
  
  v1 = round(100*v1, 0.1);
  v2 = round(100*v2, 0.1);
  v3 = round(100*v3, 0.1);
  v4 = round(100*v4, 0.1);
  v5 = round(100*v5, 0.1);
  
  
run;

proc report data = ss2 nowd split="#";
  columns rdrate groupn group typen type v1-v5;
  define rdrate / order=data "Missingness#Scenario";
  define groupn / order=data noprint;
  define group  / order=data "Treatment#Group" ;
  define typen  / order=data noprint;
  define type / order=data "Type"  ;
  define v1 / order=data "Visit 1" display format=8.1;
  define v2 / order=data "Visit 2" display format=8.1;
  define v3 / order=data "Visit 3" display format=8.1;
  define v4 / order=data "Visit 4" display format=8.1;
  define v5 / order=data "Visit 5" display format=8.1;

run;
  
data ss3;
  set ss1;

  *** PROPORTION OF OBSERVED DISCONTINUED ***;  
  rx1 = x1 / (1-o1);
  rx2 = x2 / (1-o2);
  rx3 = x3 / (1-o3);
  rx4 = x4 / (1-o4);
  rx5 = x5 / (1-o5);

  *** PROPORTION OF WITHDRAWN DISCONTINUED ***;
  rm1 = m1 / (1-o1);
  rm2 = m2 / (1-o2);
  rm3 = m3 / (1-o3);
  rm4 = m4 / (1-o4);
  rm5 = m5 / (1-o5);
  
run;  
  
  
data ss4 (keep = rdrate group: type: v:);
  set ss3 (in = in1 keep = groupn group rdrate rx1-rx5 rename = (rx1=v1 rx2=v2 rx3=v3 rx4=v4 rx5=v5))
      ss3 (in = in2 keep = groupn group rdrate rm1-rm5 rename = (rm1=v1 rm2=v2 rm3=v3 rm4=v4 rm5=v5));
  by rdrate groupn group;
  
  if in1 then typen = 1;
  else if in2 then typen = 2;
  
  select (typen);
    when (1) type = "O";
    when (2) type = "X";
    otherwise;
  end;
  
  v1 = round(100*v1, 0.1);
  v2 = round(100*v2, 0.1);
  v3 = round(100*v3, 0.1);
  v4 = round(100*v4, 0.1);
  v5 = round(100*v5, 0.1);
  
  
run;

proc report data = ss4 nowd split="#";
  columns rdrate groupn group typen type v1-v5;
  define rdrate / order=data "Missingness#Scenario";
  define groupn / order=data noprint;
  define group  / order=data "Treatment#Group" ;
  define typen  / order=data noprint;
  define type / order=data "Type"  ;
  define v1 / order=data "Visit 1" display format=8.1;
  define v2 / order=data "Visit 2" display format=8.1;
  define v3 / order=data "Visit 3" display format=8.1;
  define v4 / order=data "Visit 4" display format=8.1;
  define v5 / order=data "Visit 5" display format=8.1;

run;
  
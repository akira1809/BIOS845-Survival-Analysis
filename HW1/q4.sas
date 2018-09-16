dm 'log; clear; output; clear;';


proc format;
	value status 0 = "censored"
				 1 = "event";
	value group  1 = "AZT + ddc"
				 2 = "AZT + ddc + saquinivir";
;

data Q4;
	input group time status @@;
	format group group.
	       status status.;
	datalines;
	1 4 0 	1 6 1	1 11 1	1 12 1	1 32 1 	1 35 1
	1 38 0 	1 39 1	1 45 1	1 49 1	1 75 1	1 80 1
	1 84 1	1 85 1	1 87 1	1 102 1	1 180 0 2 2 1
	2 3 1	2 4 1	2 12 1	2 22 1	2 48 1	2 51 0
	2 56 0	2 80 1	2 85 1	2 90 1	2 94 0	2 160 1
	2 171 1	2 180 1	2 180 0	2 238 1
	;
run;

proc print data=Q4;
run;

title "Question 4 KM Estimate";

proc lifetest data=Q4;
	time time*status(0);
	by group;
run;

title "Compare H hat and H tilde";
proc lifetest data=Q4 nelson;
	time time*status(0);
	by group;
	ods output productlimitestimates=ple;
run;

data cumhazard;
	set ple;
	hatcumhaz = -log(Survival);
run;

proc print data=cumhazard;
run;

title "cuulative hazard fuction estimate comparison";
title2 "nelson-aalen estimate and product limit estimate";
proc sgpanel data=cumhazard
	(where = (censor = 0 and (group = 1 or group = 2)));
	panelby group;
	step x =time y = CumHaz/legendlabel = "N-A estimates";
	step x = time y = hatCumHaz/legendlabel = "product limit estimates";
	rowaxis label = "cumulative hazard function estimate";
run;

title "life table estimate";
proc lifetest data=Q4 method = LT width = 60 plots = (H);
	time time*status(0);
	by group;
run;

title "confidence bands";
title2 "linear transformation";
proc lifetest data=Q4 
     outsurv = results1 conftype = linear confband = all 
	 plots = survival(CB = all);
	time time*status(0);
	by group;
run;

proc print data=results1;
run;


title2 "log-log transformation";
proc lifetest data=Q4 
     outsurv = results2 conftype = loglog confband = all 
	 plots = survival(CB = all);
	time time*status(0);
	by group;
run;

proc print data=results2;
run;

title2 "arcsin-square root transformation";
proc lifetest data=Q4 
     outsurv = results3 conftype = asinsqrt confband = all 
	 plots = survival(CB = all);
	time time*status(0);
	by group;
run;

proc print data=results3;
run;


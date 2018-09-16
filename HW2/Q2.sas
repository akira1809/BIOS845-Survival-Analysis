dm 'log; clear; output; clear;';

proc format;
	value status 0 = "censored"
				 1 = "event";
	value group  1 = "AZT + ddc"
				 2 = "AZT + ddc + saquinivir";
;

data Q2;
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

/*part (a): logrank test*/
/*part (b): wilcoxon test*/

proc lifetest data=Q2;
	time time*status(0);
	strata /group = group test = (logrank fleming(0, 0.5));
run;

/*part (c): exponential AFT model*/
/*part (d): weibull AFT model*/
proc lifereg data=Q2;
	class group;
	model time*status(0) = group
	/dist = exponential;
run;

proc lifereg data=Q2;
	class group;
	model time*status(0) = group
	/dist = weibull;
run;

/*part (e) goodness of fit test*/

data fit;
	pvalue = 1 - probchi(0.117,1);
run;

proc print data=fit;
run;

data fit2;
	pvalue2 = 1 - probchi(295.533 - 295.416, 1);
run;

proc print data=fit2;
run;

proc lifereg data=Q2;
	class group;
	model time*status(0) = group
	/dist = gamma;
run;

data fit3;
	pvalue1 =  1 - probchi(98.426 - 96.584, 2);
	pvalue2 = 1 - probchi(98.309 - 96.584, 1);
run;

proc print data=fit3;
run;

/*(f) plot log survival and log log survival to check goodness of fit*/
proc lifetest data=Q2 plots(only) = ls;
	time time*status(0);
run;

proc lifetest data=Q2 plots(only) = lls;
	time time*status(0);
run;

/*(g) fit lognormal and log-logistic and report log-likelihood values*/
proc lifereg data=Q2;
	class group;
	model time*status(0) = group
	/dist = lnormal;
run;

proc lifereg data=Q2;
	class group;
	model time*status(0) = group
	/dist = llogistic;
run;

/*(i) sketch the survival curves*/
data one;
	C = 1; group = 1; output;
	C = 1; group = 2; output;
run;

data Q2_2;
	set Q2 one;
run;

proc print data=Q2_2;
run;

proc lifereg data=Q2_2;
	class group;
	model time*status(0) = group
	/dist = exponential;
	output out = surv_est quantiles = 0.02 to 0.98 by 0.02
	predicted = pred control = C;
run;

data surv_est2;
	set surv_est; sdf = 100*(1 - _prob_);
run;

proc print data=surv_est2;
run;

proc gplot data=surv_est2;
	plot sdf*pred=group;
	symbol1 I = spline color = red L= 1;
	symbol2 I = spline color = blue L = 1;
run;





dm 'log;clear;output;clear';

libname data "C:\akira\data";

proc format;
	value treatment 0 = "CPVM"
				    1 = "BCG";
run;

/*Q1 part a*/
proc lifetest data=data.Melanoma plots=s;
	time surv*dead(0);
	strata trt;
run;

/*Q1 part b*/
proc sort data=data.Melanoma out=Melanoma_sorted;
	by surv;
run;

proc sort data=data.Melanoma out=Melanoma_sorted2;
	by trt;
run;

proc print data=Melanoma_sorted;
run;

proc print data=Melanoma_sorted2;
run;

/*Q1 part c*/
proc lifetest data= data.Melanoma;
	time surv*dead(0);
	strata trt/test = (logrank wilcoxon);
run;

proc lifetest data= data.Melanoma;
	time surv*dead(0);
	strata /group = trt test = (logrank wilcoxon);
run;



/*Q1 part d*/
proc lifetest data= data.Melanoma;
	time remdur*recur(0);
	strata trt/test = (logrank wilcoxon);
run;

proc lifetest data= data.Melanoma;
	time remdur*recur(0);
	strata /group = trt test = (logrank wilcoxon);
run;

/*Q1 part E*/

proc lifetest data=data.Melanoma;
	time surv*dead(0);
	test age  sex  trt;
run;

/*Q1 part F*/
proc lifetest data=data.Melanoma plots=(ls lls);
	time surv*dead(0);
run;

/*Q2*/
/*exponential AFT*/
proc lifereg data=data.Melanoma;
	class trt; 
	model surv*dead(0) = trt sex age/dist = exponential;
run;

/*weibullAFT*/
proc lifereg data=data.Melanoma;
	class trt; 
	model surv*dead(0) = trt sex age/dist = weibull;
run;

/*generalized gamma AFT*/
proc lifereg data=data.Melanoma;
	class trt; 
	model surv*dead(0) = trt sex age/dist = gamma;
run;

/*log-normal AFT*/
proc lifereg data=data.Melanoma;
	class trt; 
	model surv*dead(0) = trt sex age/dist = lnormal;
run;

/*log-logistic AFT*/
proc lifereg data=data.Melanoma;
	class trt; 
	model surv*dead(0) = trt sex age/dist = llogistic;
run;

/*goodness of fit for nested models: exp, weibull, generalized gamma*/
data fit;
	pvalue_gamma_vs_exp =1 - probchi(46.360 - 38.283, 2);
	pvalue_gamma_vs_weib= 1 - probchi(43.474 -38.283, 1); 
	pvalue_weib_vs_expo = 1 - probchi(46.360 - 43.474, 1);
run;

proc print data=fit;
run;

title "NLMIXED: Gompertz distribution"; 
proc nlmixed data=data.Melanoma; 
	parms log_gamma -2; 
	gamma = exp(log_gamma); 
    linp = b0 + b1*trt + b2*sex + b3*age; 
	alpha = exp(-linp); 
	G_t = exp((alpha/gamma)*(1 - exp(gamma*surv))); 
	g = alpha*exp(gamma*surv)*G_t; 
	ll = (dead=1)*log(g) + /* ll for observed failures */ 
		(dead=0)*log(G_t); /* ll for censored failures */ 
	model ll ~ general(ll); 
	estimate "gamma" exp(log_gamma); 
run;

proc nlmixed data=data.Melanoma; 
	parms log_gamma -9; 
	gamma = exp(log_gamma); 
    linp = b0 + b1*trt + b2*sex + b3*age; 
	alpha = exp(-linp); 
	G_t = exp((alpha/gamma)*(1 - exp(gamma*surv))); 
	g = alpha*exp(gamma*surv)*G_t; 
	ll = (dead=1)*log(g) + /* ll for observed failures */ 
		(dead=0)*log(G_t); /* ll for censored failures */ 
	model ll ~ general(ll); 
	estimate "gamma" exp(log_gamma); 
run;

data Melanoma;
	set data.melanoma;
	format trt treatment.;
run;

proc print data= melanoma;
run;

proc lifereg data=Melanoma;
	class trt; 
	model surv*dead(0) = trt sex age/dist = lnormal;
run;

proc lifereg data=data.Melanoma;
	class trt; 
	model surv*dead(0) = trt sex age/dist = lnormal;
run;

/*Q2 part B*/
title "survival curve for female age 40 by treatment group";
data one;
	C = 1; trt = 0; age = 40; sex = 1; output;
	C = 1; trt = 1; age = 40; sex = 1; output;
run;

data Melonoma;
	set data.Melanoma one;
run;

proc print data=Melonoma;
run;

proc lifereg data=Melonoma;
	class trt;
	model surv*dead(0) = trt sex age
	/dist = lognormal;
	output out = surv_est quantiles = 0.02 to 0.98 by 0.02
	predicted = pred control = C;
run;

proc print data=surv_est;
run;

data surv_est2;
	set surv_est; sdf = 100*(1 - _prob_);
run;

proc print data=surv_est2;
run;

proc gplot data=surv_est2;
	plot sdf*pred=trt;
	symbol1 I = spline color = red L= 1;
	symbol2 I = spline color = blue L = 1;
run;

/*Q3*/

data Addicts;
	set data.Addicts;
	where status = 1;
run;

proc freq data=Addicts noprint;
	table length
	/out = ties;
run;

proc print data=ties;
	 where count >= 2;
run;

proc freq data=data.Addicts;
	table status;
run;

/*part C, default tie handling with breslow*/
proc phreg data=data.Addicts;
	model length*status(0) = clinic prison dose/risklimits;
run;

/*part D other tie handling*/
proc phreg data=data.Addicts;
	model length*status(0) = clinic prison dose/ties = exact risklimits;
run;

proc phreg data=data.Addicts;
	model length*status(0) = clinic prison dose/ties = efron risklimits;
run;

/*part E make dosage a categorical variable*/
data Addicts;
	set data.Addicts;
	if (dose < 60) then do;
		dose1 = 0;
		dose2 =0;
	end;
	else if (60 <= dose < 80) then do;
		dose1 = 1;
		dose2 = 0;
	end;
	else if (dose >=80) then do;
		dose1 = 0;
		dose2 = 1;
	end;
run;

proc print data=Addicts;
run;	

proc phreg data=Addicts;
	class dose1(ref = "0") dose2(ref = "0") clinic(ref = "0");
	model length*status(0) = clinic prison dose1 dose2 /risklimits;
run;

/*part F*/
proc lifetest data=Addicts plots= (ls lls);
	time length*status(0);
run;

/*part G*/
proc lifetest data=Addicts plots = lls;
	time length*status(0);
	strata clinic/;
run;

proc lifetest data=Addicts plots = lls;
	time length*status(0);
	strata prison/;
run;

proc sort data=addicts out = addicts_sorted;
	by length;
run;

proc print data=addicts_sorted;
run;

proc lifetest data=Addicts method = life outsurv = results
		plots = (H) intervals = 50 100 200 300 400 500 600 700 800 900 1000 1100;
	time length*status(0);
	strata clinic;
run;

data loghazard;
	set results;
	loghazard = log(hazard);
run;

proc sgplot data=loghazard;
	series x = length y = loghazard/group = clinic;
run;

proc lifetest data=Addicts method = life outsurv = results2
		plots = (H) intervals = 50 100 200 300 400 500 600 700 800 900 1000 1100;
	time length*status(0);
	strata prison;
run;

data loghazard2;
	set results2;
	loghazard2 = log(hazard);
run;

proc sgplot data=loghazard2;
	series x = length y = loghazard2/group = prison;
run;


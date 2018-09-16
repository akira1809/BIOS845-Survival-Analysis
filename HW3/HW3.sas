dm 'log; clear; output; clear;';

libname data "C:\akira\data";

/*Question 1*/
/*part A fit a cox model with age, sex, age*sex*/

data Whas_1;
	set data.Whas;
	agesex = age * sex;
run;

proc phreg data=Whas_1;
	model LENFOL*FSTAT(0) = age sex agesex;
run;

proc phreg data=Whas_1;
 	model LENFOL*FSTAT(0) = age sex agesex/selection = FORWARD;
run;

proc phreg data=Whas_1;
 	model LENFOL*FSTAT(0) = age sex agesex/selection = BACKWARD;
run;

proc phreg data=Whas_1;
 	model LENFOL*FSTAT(0) = age sex agesex/selection = STEPWISE;
run;

/*part C compute and graph estimated survival functions 
for 65 year old males and females, and estimate median survival times*/

data groups;
	input age sex agesex;
	datalines;
	65 0 0
	65 1 65
	;
run;

proc phreg data=Whas_1;
	model LENFOL*FSTAT(0) = age sex agesex;
	baseline covariates = groups out=cov_adj_surv
		survival = surv;
run;

proc print data=cov_adj_surv noobs;
	by sex;
run;

proc gplot data=cov_adj_surv;
	plot surv*LENFOL = sex;
	symbol1 I = STEPLJ COLOR = BLUE L = 1;
	symbol2 I = STEPLJ COLOR = Red L = 2;
run;

/*part D fit a cox model with agecat, sex, agecat*sex*/
data Whas_2;
	set data.Whas;
	IF age_cat = 1 THEN cat1 = 1;
	ELSE cat1 = 0;
	IF age_cat = 2 THEN cat2 = 1;
	ELSE cat2 = 0;
	IF age_cat = 3 THEN cat3 = 1;
	ELSE cat3 = 0;
	cat1_sex = cat1*sex;
	cat2_sex = cat2*sex;
	cat3_sex = cat3*sex;
run;

proc print data=Whas_2;
run;

proc phreg data=Whas_2;
	model LENFOL*FSTAT(0) = cat1 cat2 cat3 sex cat1_sex cat2_sex cat3_sex;
run;

proc phreg data=Whas_2;
	model LENFOL*FSTAT(0) = cat1 cat2 cat3 sex cat1_sex cat2_sex cat3_sex/selection = stepwise;
run;

proc phreg data=Whas_2;
	model LENFOL*FSTAT(0) = cat1 cat2 cat3 sex cat1_sex cat2_sex cat3_sex/selection = forward;
run;

proc phreg data=Whas_2;
	model LENFOL*FSTAT(0) = cat1 cat2 cat3 sex cat1_sex cat2_sex cat3_sex/selection = backward;
run;


proc phreg data=data.Whas;
	class age_cat;
	model LENFOL*FSTAT(0) = age_cat sex age_cat*sex;
run;

/*part E*/
proc phreg data=Whas_2;
	model LENFOL*FSTAT(0) = cat1 cat2 cat3 sex cat1_sex cat2_sex cat3_sex/risklimits;
run;

/*part F validate PH model*/
/*model include sex and age*/

/*use Schoenfeld Residual test checking on age*/
/*no apparent pattern, PH model valid for age*/
proc phreg data=data.Whas;
	model LENFOL*FSTAT(0) = age sex;
	output out=resid1
		RESSCH = SCHAGE SCHSEX;
run;

proc gplot data=resid1;
	plot SCHAGE*LENFOL;
	symbol1 value = circle H=.5;
run;

/*use log-log survival functions to check on sex*/
/*two curves parallel, PH model valid for sex*/
proc lifetest data=data.Whas notable
	plots = (LLS) graphics cs=none;
	time LENFOL*FSTAT(0);
	strata sex;
run;


/*use observed versus expected plots of survival functions
to check on sex*/
/*for both sex, the observed and expected plots are very close,
PH model is valid for sex.*/

proc lifetest data=data.Whas notable
	outsurv = observed;
	time LENFOL*FSTAT(0);
	strata sex;
run;

data covs;
	input sex @@;
	datalines;
	0 1
	;
run;

proc phreg data=data.Whas;
	model LENFOL*FSTAT(0) = sex;
	baseline covariates = COVS out=expected
		survival = survival/nomean;
run;


data combined;
	set observed(IN = IN1) expected(IN = IN2);
	if IN1 Then grp = sex;
	else if IN2 then grp = sex + 2;
run;

proc print data=combined;
run;

proc gplot data=combined;
	plot survival*LENFOL = grp;
	symbol1 I = STEPLJ COLOR = RED L = 1;
	symbol2 I = STEPLJ COLOR = BLUE L = 1;
	symbol3 I = STEPLJ COLOR = RED L = 2 W = 2;
	symbol4 I = STEPLJ COLOR = BLUE L = 2 W = 2;
run;

/*use interaction with time to check on age and sex
we add linear interaction term here*/
/*both interaction terms are not significant, appear to be no
violation of PH assumption*/

proc phreg data=data.Whas;
	model LENFOL*FSTAT(0) = age sex aget sext;
	aget = age*LENFOL;
	sext = sex*LENFOL;
run;


/*part G testing goodness of fit, outliers,
influential observations*/

/*log survival plot of coxsnell residuals*/
proc phreg data=data.Whas;
	model LENFOL*FSTAT(0) = age sex;
	output out=resid logsurv = coxsnell;
run;

data resid2;
	set resid;
	coxsnell=-coxsnell;
run;

proc lifetest data=resid2 plots = (LS) notable;
	time coxsnell*FSTAT(0);
run;

/*use martingale residual to check on functional form of age*/
proc phreg data=data.Whas;
	class sex;
	model LENFOL*FSTAT(0) = age sex;
	output out=MARTIN RESMART = MARTINGALE;
run;

/*no particular pattern. Current functional form of age seems 
to be appropriate*/
proc sgplot data=MARTIN;
	title "MARTINGALE RESIDUALS vs AGE";
	LOESS Y = MARTINGALE X = age;
run;

/*observed cumulative martingale residuals within the spectrum
of the simulated realizations, do not see problem using age
in its current form*/
proc phreg data=data.Whas;
	class sex;
	model LENFOL*FSTAT(0) = age sex;
	ASSESS var = (age)/CRPANEL;
run;

/*use deviance residuals to identify outliers*/
proc phreg data=data.Whas;
	model LENFOL*FSTAT(0) = age sex;
	output out=RESID3 RESDEV = DEV;
run;

proc univariate data=RESID3;
	var dev;
run;

title "Deviance Residuals vs Age";
proc gplot data=RESID3;
	plot DEV*AGE;
	SYMBOL1 VALUE=circle H=.5;
run;

/*use DFBETA keyword to find influential observations*/
proc phreg data=data.Whas;
	model LENFOL*FSTAT(0) = age sex;
	output out=RESID4 DFBETA=beta_age beta_sex;
run;

data resid4;
	set resid4;
	obs = _n_;
run;

proc print data=resid4;
run;

proc gplot data=resid4;
	title "influence of beta-age vs observation number";
	plot beta_age*obs;
	symbol1 value=dot H = .5;
run;

/*testing*/
/*proc phreg data=data.leukemiab;*/
/*	model weeks*relapse(0) = rx logwbc sex;*/
/*	output out=RESID5 DFBETA=beta_rx beta_lwbc beta_sex;*/
/*run;*/
/**/
/*data resid5;*/
/*	set resid5;*/
/*	obs = _n_;*/
/*run;*/
/**/
/*proc print data=resid5;*/
/*run;*/

/*Question 2*/
/*Analysis on Pneumonia data*/

/*part A treating drug as continuous covariate*/
/*adjusting for age*/

title "pneumonia cure time analysis";
proc phreg data = data.Pneumonia;
	model time*cured(0) = age drug;
run;

/*part B treating a new variable recording the 
remaining drug dosage by end of study: 
remain = initial dosage *exp(-0.35t)*/

data Pneumonia_b;
	set data.Pneumonia;
	remain = drug*exp(-0.35*time);
run;

proc print data=Pneumonia_b;
run;

proc phreg data=Pneumonia_b;
	model time*cured(0) = age remain;
run;

/*part C consider different effects over time for 
different initial dossages*/

proc phreg data=data.Pneumonia;
	model time*cured(0) = age remain;
	remain = 0;
	remain = drug * exp(-0.35*time);
run;


/*Question 3*/

/*part C fit the model developed in part A*/
/*give estimate of the hazard ratios in part B*/
proc phreg data=data.Leukemiab;
	model weeks*relapse(0)= RX LOGWBC SEX*TIME;
	TIME = 0;
	IF weeks < 4 then TIME= 1;
	ELSE IF 4 <= weeks < 8  then TIME = 3;
	ELSE IF 8 <= weeks < 12 then TIME = 5;
	ELSE IF 12 <= weeks < 16 then TIME = 7;
	ELSE IF weeks >= 16 then TIME = 9; 
run;

/*test nyh's code*/
/*proc phreg data=data.Leukemiab;*/
/*	model weeks*relapse(0)= RX LOGWBC sext;*/
/*	IF (weeks < 4) then sext= sex*1;*/
/*	ELSE IF (4 <= weeks < 8)  then sext = sex*3;*/
/*	ELSE IF (8 <= weeks < 12) then sext = sex*5;*/
/*	ELSE IF (12 <= weeks < 16) then sext = sex*7;*/
/*	ELSE IF (weeks >= 16) then sext = sex*9; */
/*run;*/

/*part D fit a cox model stratified on sex*/
proc sort data=data.Leukemiab out=Leukemiab_sorted;
	by SEX;
run;

proc phreg data=Leukemiab_sorted;
	by SEX;
	model weeks*relapse(0)= RX LOGWBC;
run;

/*Question 5*/

/*part B*/
/*(i) run the code given on page 165*/
data RECIDCUM;
	set data.Recid;
	ARRAY emp(*) emp1-emp52;
	ARRAY cum(*) cum1-cum52;
	cum1 = emp1;
	DO i = 2 to 52;
		cum(i) = (cum(i - 1)*(i - 1) + emp(i))/i;
	END;
run;

proc phreg data=recidcum;
	where week > 1;
	model week*arrest(0) = fin age race wexp mar paro prio employ/TIEs = EFRON;
	array cumemp(*) cum1-cum52;
	employ = cumemp[week-1];
run;

/*(ii) code inside proc phreg to do the same 
job as in (i)*/

proc phreg data=data.Recid;
	where week > 1;
	model week*arrest(0) = fin age race wexp mar paro prio employ/TIEs = EFRON;
	array emp[*] emp1-emp52;
	array cum[*] cum1-cum52;
	cum1 = emp1;
	DO i = 2 to 52;
		cum[i] = (cum[i - 1]*(i - 1) + emp[i])/i;
		IF week = i THEN employ = cum[i-1];
	END;
run;	

/*part C*/
/*(i) consider number of switches as a time dependent 
covariate*/
proc phreg data=data.Recid;
	model week*arrest(0) = fin age race wexp mar paro prio num_switch/TIEs = EFRON;
	array emp[*] emp1-emp52;
	array switch[*] switch1-switch52;
	switch1 = 0;
	num_switch = 0;
	Do i = 2 to 52;
		IF emp[i] - emp[i -1] ne 0 THEN switch[i] = 1;
		ELSE switch[i] = 0;
		num_switch = num_switch + switch[i];
	END;
run;

/*(ii) consider number of negative switches as a time 
dependent covariate*/
proc phreg data=data.Recid;
	model week*arrest(0) = fin age race wexp mar paro prio num_switch/TIEs = EFRON;
	array emp[*] emp1-emp52;
	array switch[*] switch1-switch52;
	switch1 = 0;
	num_switch = 0;
	Do i = 2 to 52;
		IF emp[i] - emp[i -1] =-1 THEN switch[i] = 1;
		ELSE switch[i] = 0;
		num_switch = num_switch + switch[i];
	END;
run;

/*Question 6*/

/*continue from 5(C)(ii), categorize number of negative switch
as a time dependent covariate*/

/*switchtype_1 is a time dependent dummy covariate indicating #of job lost = 0*/
/*switchtype_2 is a time dependent dummy covariate indicating #of job lost = 1*/

proc phreg data=data.Recid;
	model week*arrest(0) = fin age race wexp mar paro prio switchtype_1 switchtype_2/TIEs = EFRON;
	array emp[*] emp1-emp52;
	array switch[*] switch1-switch52;
	switch1 = 0;
	num_switch = 0;
	Do i = 2 to 52;
		IF emp[i] - emp[i -1] =-1 THEN switch[i] = 1;
		ELSE switch[i] = 0;
		num_switch = num_switch + switch[i];
		IF num_switch = 0 THEN switchtype_1 = 1;
		ELSE switchtype_1 = 0;
		IF num_switch = 1 THEN switchtype_2 = 1;
		ELSE switchtype_2 = 0;
	END;
run;


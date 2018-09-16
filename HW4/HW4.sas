dm 'log; clear; output; clear;';

libname data "C:\akira\data";

/*HW4*/

/*Question 1*/

/*part A*/
/*account for different proportional hazard for age
before and after week 15*/

proc phreg data=data.leukemiab;
	model weeks*relapse(0) = RX LOGWBC SEX SEXT/TIES = exact risklimits;
	SEXT = SEX*(weeks > 15);
 	BEFORE15: TEST SEX =0;
	AFTER15: TEST SEX+SEXT = 0; 
run;

/*part B*/
/*fitting a complementary log-log proportional hazard model*/
/*fitting a linear effect of time, since we have empty intervals
by listing the variable on the model statement alone, instead of using 
a class statement*/
/*also consider a time dependent variable for sex before and after 15 weeks*/
data leukemiab_extend;
	set data.leukemiab;
	do period = 1 to weeks;
		IF period = weeks and relapse = 1 then response = 1;
		else response = 0;
		sext = sex*(period > 15);
		output;
	end;
run;

proc logistic data=leukemiab_extend;
	model response(event = '1') = RX LOGWBC SEX SEXT period/link = cloglog CLPARM=BOTH;
run;

/*Question 1 part II*/
/*part B*/
/*fit the complementary log-log model to the job duration data*/
data jobyrs;
	set data.jobdur;
	do year = 1 to dur;
		IF year = dur and event = 1 then response = 1;
		else response = 0;
		log_year = log(year);
		output;
	end;
run;

proc logistic data=jobyrs;
	class year/param = glm;
	model response(event = '1') = ed prestige salary year/link = cloglog;
run;

proc logistic data=jobyrs;
	model response(event = '1') = ed prestige salary year/link = cloglog;
run;

proc logistic data=jobyrs;
	model response(event = '1') = ed prestige salary year year*year/link = cloglog;
run;

proc logistic data=jobyrs;
	model response(event = '1') = ed prestige salary log_year/link = cloglog;
run;

/*Question 2*/
/*part I*/

/*part A*/
/*we assume also 3 years accrual time
and 1 year follow up, then it is doable with SAS*/
proc power;
	twosamplesurvival test = logrank
	curve("control") = (1):(0.93)
	refsurvival = "control"
	hazardratio = 1.38
	accrualtime = 3
	followuptime = 1
	eventstotal = .
	power= 0.9
	;
run;

/*part B*/
proc power;
	twosamplesurvival test = logrank
	curve("control") = (1):(0.93)
	refsurvival = "control"
	hazardratio = 1.38
	accrualtime = 3
	followuptime = 1
	ntotal = .
	power= 0.9
	;
run;
/*equivalent as the code above, only that here
we compute sample size in each group*/
proc power;
	twosamplesurvival test = logrank
	curve("control") = (1):(0.93)
	refsurvival = "control"
	hazardratio = 1.38
	accrualtime = 3
	followuptime = 1
	npergroup = .
	power= 0.9
	;
run;

/*part C*/
proc power;
	twosamplesurvival test = logrank
	curve("control") = (1):(0.93)
	refsurvival = "control"
	hazardratio = 1.38
	accrualtime = 3
	followuptime = 1
	power= .
	ntotal = 2000;
	plot x = n min =500 max = 3000;
run;

/*part D*/
/*group weight: control vs treatment = 1:3*/
/*events total*/
proc power;
	twosamplesurvival test = logrank
	curve("control") = (1):(0.93)
	refsurvival = "control"
	hazardratio = 1.38
	accrualtime = 3
	followuptime = 1
	groupweights = (1 3)
	eventstotal = .
	power= 0.9
	;
run;

/*total sample size*/
proc power;
	twosamplesurvival test = logrank
	curve("control") = (1):(0.93)
	refsurvival = "control"
	hazardratio = 1.38
	accrualtime = 3
	followuptime = 1
	groupweights = (1 3)
	ntotal = .
	power= 0.9
	;
run;


/*part II of Question 2*/
proc power;
	twosamplesurvival test = logrank
	curve("control") = (0.25 0.5 0.75 1 1.2 1.5 2 2.75 6)
						:(0.95 0.9 0.65 0.5 0.375 0.225 0.15 0.075 0.05)
	curve("treatment") = (0.5 1 1.5 2 2.4 3 4 5.5 12)
						:(0.95 0.9 0.65 0.5 0.375 0.225 0.15 0.075 0.05)
	gsurv = "control"|"treatment"
	accrualtime = 3
	followuptime = 2
	npergroup = .
	power = 0.8
	alpha = 0.1
	sides = 1;
run;

/*Question 3*/
/*Part I Exercise 5.2*/

/*we follow the algorithm step by step*/

dm 'log; clear; output; clear;';
proc iml;
	/*input our data and set up initial values for important vectors*/
	/*t is the time grid*/
	t = {0, 2, 4, 8, 10, 12, 15, 20, 24, 28, 34, 41, 62, 69, 75, 79, 86};
	/*initialize survival. S[1] is for t = 0 so it is always 1*/
	S = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	m = 17;
	/*c is the left censoring number at each time*/
	c = {0, 0, 0, 0, 1, 2, 2, 4, 3, 3, 4, 2, 3, 2, 1, 2, 3};
	/*d is the number of events at each time*/
	d = {0, 3, 2, 1, 2, 4, 6, 3, 3, 2, 1, 0, 0, 0, 0, 0, 0};
	/*r is the number of right censoring at each time*/
	r = {0, 0, 0, 0, 0, 0, 1, 1, 2, 3, 5, 3, 4, 6, 6, 3, 7};

	/*step 0: produce an initial estimate of the survival function S_0(t_j)
by using product limits estimate (KM)*/
	/*computing risking set, ignore left censoring*/
	Y = {68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	do i = 2 to m;
		do j = i to m;
			Y[i] = Y[i] + d[j] + r[j];
		end;
	end;
	print y;
	do j = 2 to m;
		S[j] = S[j - 1]*(Y[j] - d[j])/Y[j];
	end;
	/*we checked that we have the same result as from proc lifetest for step 0*/
	print s;

/*k is a pseudo index here. just to make the same program run 3 times*/
do k = 1 to 3;
	/*Step (k)1- (k)2: estimate p_ij, estimate the number of events at time t_j*/
	do j = 2 to m;
		do i = j to m;
			d[j] = round(d[j] + c[i]*(s[j - 1] - s[j])/(1 - s[i]), .0001);
		end;
	end;
	/*Step (k)3 compute the product limit estimate as in step 0, using the new d*/
	Y = {68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	do i = 2 to m;
		do j = i to m;
			Y[i] = Y[i] + d[j] + r[j];
		end;
	end;
	do j = 2 to m;
		S[j] = round(S[j - 1]*(Y[j] - d[j])/Y[j], .0001);
	end;
end;
print s;

quit;

/*test*/
/*data test;*/
/*	input time status @@;*/
/*	datalines;*/
/*	2 1 	2 1 	2 1		4 1		4 1		8 1		10 1	10 1*/
/*	12 1	12 1	12 1	12 1	15 1	15 1	15 1	15 1*/
/*	15 1	15 1	15 0	20 1	20 1	20 1	20 0	24 1*/
/*	24 1	24 1	24 0	24 0	28 1	28 1	28 0	28 0*/
/*	28 0	34 1	34 0	34 0	34 0	34 0	34 0	41 0*/
/*	41 0	41 0	62 0	62 0	62 0	62 0	69 0	69 0*/
/*	69 0	69 0	69 0	69 0	75 0	75 0	75 0	75 0*/
/*	75 0	75 0	79 0	79 0	79 0	86 0	86 0	86 0*/
/*	86 0	86 0	86 0	86 0*/
/*	;*/
/*run;*/
/**/
/*proc lifetest data=test outsurv = test_result noprint;*/
/*	time time*status(0);*/
/*run;*/
/**/
/*proc print data=test_result;*/
/*	var time survival;*/
/*run;*/

/*Question 3 part II*/

/*reformat the data*/
data ex5_2;
	input ltime rtime @@;
	datalines;
	2 2 	2 2 	2 2 	4 4 	4 4		8 8 
	10 10 	10 10 	. 10	12 12 	12 12 	12 12
	12 12	. 12	. 12	15 15 	15 15 	15 15
	15 15	15 15	15 15	. 15	. 15	15 .
	20 20	20 20	20 20	. 20	. 20	. 20
	. 20	20 .	24 24 	24 24 	24 24	. 24
	. 24	. 24	24 .	24 .	28 28	28 28
	. 28	. 28	. 28	28 .	28 .	28 .
	34 34 	. 34	. 34	. 34	. 34	34 .
	34 .	34 .	34 .	34 .	. 41	. 41
	41 .	41 .	41 .	. 62	. 62	. 62
	62 .	62 .	62 .	62 .	. 69	. 69
	69 .	69 .	69 .	69 .	69 .	69 .
	. 75	75 .	75 .	75 .	75 .	75 .
	75 .	. 79	. 79	79 .	79 .	79 .
	. 86	. 86	. 86	86 .	86 .	86 .
	86 .	86 .	86 .	86 .
	;
run;

proc iclifetest data=ex5_2 impute(seed = 1234) outsurv = compare;
	time (ltime rtime);
run;

proc print data=compare;
	var leftboundary rightboundary survprob;
run;


/*Question 4*/
/*input the data*/

data Q4;
	/*for group, 1 = adult, 0 = children	*/
	input group left right @@;
	datalines;
	1 0 4	1 0 24	1 0 24	1 0 36	1 0 36	1 0 36	1 4 8	1 4 8
	1 4 8	1 4 8	1 4 8	1 4 24	1 4 24	1 4 24	1 4 36	1 4 48
	1 8 12	1 8 12	1 8 12	1 8 12	1 8 36	1 12 24 1 12 24	1 12 24
	1 12 24	1 12 24	1 12 48	1 12 48 1 12 48	1 12 48	1 24 36	1 24 36
	1 24 36	1 24 36	1 36 48	1 36 48	1 36 48	1 36 48	1 48 .	1 48 .
	1 48 .	1 48 .	1 48 .	1 48 .	1 48 .	1 48 .	 0 0 12	 0 0 12
	0 8 12	0 8 36	0 8 36	0 8 48	0 12 24	0 12 48	0 24 36	0 24 36
	0 24 36	0 24 36	0 24 36	0 24 36	0 36 48	0 36 48	0 36 48	0 36 48
	0 36 48 0 36 48	0 36 48	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .
	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .
	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .
	0 48 .	0 48 .	0 48 .	0 48 .	0 48 .
	;
run;

data cov;
	group = 0; output;
	group = 1; output;
run;


/*fit the piecewise exponential model*/
proc icphreg data= Q4;
	class group;
	model (left, right) = group/basehaz = piecewiseexponential;
	hazardratio group;
	baseline covariates = cov;
run;

/*fit the cubic spline model*/

proc icphreg data=Q4;
	class group;
	model (left, right) = group/basehaz = splines;
	hazardratio group;
	baseline covariates = cov;
run;



	
	
	

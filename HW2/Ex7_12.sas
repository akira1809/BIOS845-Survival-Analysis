dm 'log; clear; output; clear;';

* Macro: _RENYI *;
* Author: Mat Davis *;
* Affiliation: University of Pennsylvania *;
* Email: davismat@mail.med.upenn.edu *;
* Description: This macro will obtain the Renyi Test Statistic *;
* testing the null hypothesis that the hazard rates *;
* for two stratified groups of subjects are equal. *;
*************************************************************************;
* MACRO PARAMETERS *;
* *;
* DATA: The dataset to be analyzed *;
* STRATA: The variable that indicates the group to be stratified *;
* TIME: The variable that indicates survival time *;
* CENSOR: The variable that indicates censored status *;
* CVAL: The value of the censored observations in CENSOR *;
* WEIGHT: The desired weight for the analysis *;
* LOGRANK: 1 *;
* GEHAN: Combined at risk set *;
* TARONEWARE: The square root of the *;
* combined risk set at time t *;
* PETOPETO: The combined survival rate at *;
* time t *;
* MODPETOPETO: Corrected survival rate at *;
* time t *;
* FH: (s(t*))^p*(1-s(t*))^q where s(t) is the *;
* combined survival rate at time t-1 *;
* FH1: The value of 'p' for the Flemming-Harrington weights *;
* FH2: The value of 'q' for the Flemming-Harrington weights *;
*************************************************************************;
%macro _renyi(data,strata,time,censor,cval,weight,fh1,fh2);
%macro _fillgap(old,new);
if not missing(&old.) then &new.=&old.;
%mend;
proc lifetest method=km plots=(s) data=&data.;
strata &strata.;
time &time.*&censor.(&cval.);
ods output productlimitestimates=prod;
run;
proc means data=prod noprint;
class &time. &strata.;
var failed;
output out=prodmax(where=(not missing(&time.) and not missing(&strata.))) max=failed;
run;quit;
proc means data=prod noprint;
class &time. &strata.;
var censor;
output out=cenmax(where=(not missing(&time.) and not missing(&strata.))) max=censor;
run;quit;
proc means data=prod noprint;
class &time. &strata.;
var left;
output out=prodmin(where=(not missing(&time.) and not missing(&strata.))) min=left;
run;quit;
data prod;
merge prodmax(keep=&time. &strata. failed) prodmin(keep=&time. &strata. left) cenmax(keep=&time. &strata. censor);
by &time. &strata.;
run;
data prod1 prod2;
set prod;
run;
data prodx(keep=&time. d1 d2 y1 y2 totf renyi1 censor1 censor2 renyisum
absrenyisum renyisd renyisdsum absrenyisdsum _id atrisk atrisk1 atrisk2
lastd1 lastd2 lastcensor1 lastcensor2 survival weight lastsurv);
retain g1 g2 f1 f2 l1 l2 renyisum 0 renyisdsum 0 atrisk atrisk1 atrisk2 survival 1;
merge prod1(where=(&strata.1=1) rename=(&strata.=&strata.1 failed=failed1 left=left1 censor=censor1))
prod2(where=(&strata.2=2) rename=(&strata.=&strata.2 failed=failed2 left=left2 censor=censor2));
by &time.;
***DATA MANIPULATION TO OBTAIN FULL TABLE***;
%_fillgap(&strata.1,g1) %_fillgap(&strata.2,g2) %_fillgap(failed1,f1)
%_fillgap(failed2,f2) %_fillgap(left1,l1) %_fillgap(left2,l2)
if missing(censor1) then censor1=0; if missing(censor2) then censor2=0;
lastf1=lag(f1); lastf2=lag(f2); last11=lag(l1); lastl1=lag(l2);
lastcensor1=lag(censor1); lastcensor2=lag(censor2);
d1=f1-lastf1; d2=f2-lastf2;
y1=l1+d1+censor1; y2=l2+d2+censor2;
lastd1=lag(d1); lastd2=lag(d2); totf=sum(d1,d2);
***REMOVE ENTRY AT TIME=0***;
if &time.=0 then delete;
***GET THE AT RISK SET***;
if _N_=2 then atrisk=sum(y1,y2);
else atrisk=atrisk-lastd1-lastd2-lastcensor1-lastcensor2;
if _N_=2 then atrisk1=y1;
else atrisk1=atrisk1-lastd1-lastcensor1;
if _N_=2 then atrisk2=y2;
else atrisk2=atrisk2-lastd2-lastcensor2;
***DETERMINE SURVIVAL RATES***;
survpart=(1-(totf/(atrisk+1)));
if not missing(survpart) then survival=survival*survpart;
lastsurv=survival/survpart;
***Correct for 0^0***;
if lastsurv=0 then lastsurv=.0000001;
if lastsurv=1 then lastsurv=.9999999;
***SELECTION OF THE WEIGHT***;
%if %upcase(&weight.)=LOGRANK %then %do;
weight=1;
%end;
%else %if %upcase(&weight.)=GEHAN %then %do;
weight=atrisk;
%end;
%else %if %upcase(&weight.)=TARONEWARE %then %do;
weight=atrisk**.5;
%end;
%else %if %upcase(&weight.)=PETOPETO %then %do;
weight=survival;
%end;
%else %if %upcase(&weight.)=MODPETOPETO %then %do;
weight=survival*atrisk/(atrisk+1);;
%end;
%else %if %upcase(&weight.)=FH %then %do;
weight=lastsurv**(&fh1.)*(1-lastsurv)**(&fh2.);
%end;
***OBTAIN THE RENYI STATISTICS***;
if atrisk=0 then renyi1=0; else renyi1=weight*(d1-atrisk1*(totf/atrisk));
renyisum=renyisum+renyi1;
absrenyisum=abs(renyisum);
if atrisk-1=0 then renyisd=0; else renyisd=weight**2*(atrisk1/atrisk)*(atrisk2/atrisk)*((atrisk-totf)/(atrisk-1))*totf;
renyisdsum=renyisd+renyisdsum;
absrenyisdsum=abs(renyisdsum);
_id=1;
run;quit;
***OBTAIN THE SUM OF THE RENYI STATISTICS (TEST STATISTIC)***;
proc means data=prodx sum noprint;
var renyi1;
output out=lr sum=lr;
run;quit;
data _null_;
set lr;
call symput('lr',lr);
run;
proc sort data=prodx;
by &time.;
run;
***GET RENYI SD***;
data _null_;
set prodx;
by _id;
if last._id;
call symput('renyisdsum',renyisdsum);
run;
proc sort data=prodx;
by descending absrenyisum &time.;
run;
***REPORT RENYI STATISTICS***;
data renyi;
set prodx;
by descending absrenyisum &time.;
if _N_=1;
renyisdeval=sqrt(&renyisdsum.);
renyistat=absrenyisum/sqrt(&renyisdsum.);
file print;
put "****************************************************************";
put "**********************RENYI STATISTICS**************************";
put "****************************************************************";
put " ";
put "Time that maximizes Z(t): " &time.;
put "%upcase(&weight.) statistic (numerator): &lr.";
put "Z(t) Statistic: " absrenyisum;
put "Z(t) Standard Deviation: " renyisdeval ;
put "Renyi Test Statistic: " renyistat;
***Initialize Infsum***;
infsum=0;
%do k=0 %to 10000;
infsum=infsum+(-1)**%eval(&k.)/(2*%eval(&k.)+1)*exp(-1*constant('PI')**2*(2*%eval(&k.)+1)**2/(8*renyistat**2));
%end;
pvalue=1-(4/constant('PI'))*infsum;
put " Approximate P-value: " pvalue;
put " ";
put "****************************************************************";
run;
%mend;

proc format;
	value status 0 = "censored"
				 1 = "event";
	value group  1 = "Aneuploid"
				 2 = "Diploid";
;

data Ex7_4;
	input group time status @@;
	format group group.
	       status status.;
	datalines;
	1 1 1	1 3 1	1 3 1	1 4 1	1 10 1	1 13 1	1 13 1	1 16 1	1 16 1
	1 24 1	1 26 1	1 27 1	1 28 1	1 30 1	1 30 1	1 32 1	1 41 1 	1 51 1
	1 65 1	1 67 1	1 70 1	1 72 1	1 73 1	1 77 1	1 91 1	1 93 1	1 96 1
	1 100 1	1 104 1	1 157 1	1 167 1	1 61 0	1 74 0	1 79 0 	1 80 0 	1 81 0
	1 87 0	1 87 0	1 88 0	1 89 0	1 93 0	1 97 0 	1 101 0	1 104 0	1 108 0
	1 109 0	1 120 0	1 131 0	1 150 0	1 231 0	1 240 0	1 400 0 2 1 1	2 3 1
	2 4 1	2 5 1	2 5 1	2 8 1	2 12 1	2 13 1	2 18 1	2 23 1	2 26 1
	2 27 1	2 30 1	2 42 1	2 56 1	2 62 1	2 69 1	2 104 1	2 104 1	2 112 1
	2 129 1	2 181 1	2 8 0	2 67 0	2 76 0	2 104 0	2 176 0	2 231 0
	;
run;

data sec1_6;
	infile "C:\akira\data\burn.txt";
	input obs Treatment Gender Race Percentage Head
	      buttock trunk upper_leg lower_leg resp_tract type
		  time_to_excision excision_ind time_to_prophy
		  prophy_ind time_to_straphy straphy_ind;
	tre = treatment + 1;
run;

/*7.12a, use Renyi statistic to run log rank test*/
%_renyi(sec1_6, tre, time_to_straphy, straphy_ind, 0, LOGRANK, 0, 0);

%_renyi(Ex7_4, group, time, status, 0, Gehan, 0, 0);


PROC IMPORT OUT= WORK.SEC1_6 
            DATAFILE= "C:\akira\school\KUMC\2018 Spring\BIOS845\HW\HW2\b
urn.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

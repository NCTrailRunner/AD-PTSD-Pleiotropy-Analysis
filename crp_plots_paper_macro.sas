

%macro mark(trta,trtb);

libname dsd "W:\Final\RESULTS\COG_UNC\PLEIO_&trta.";
options nodate nonumber;
proc append base=dsd.&trtb._fold_all data=dsd.&trtb.2p0;run;
proc append base=dsd.&trtb._fold_all data=dsd.&trtb.2p5;run;
proc append base=dsd.&trtb._fold_all data=dsd.&trtb.3p0;run;
proc append base=dsd.&trtb._fold_all data=dsd.&trtb.3p5;run;
proc append base=dsd.&trtb._fold_all data=dsd.&trtb.4p0;run;


data dsd.&trtb._fold_subs;
	set dsd.&trtb._fold_all;
	c = cat("-log10(cogaddp) >= ",group);
	if fold > 10 then delete;
	/* if lp_TICS < 6.0 then output; */
run;

proc sort data=dsd.&trtb._fold_subs;
	by group lp_cog;
run;

ods html close;
ods printer printer=tiff file="W:\Final\RESULTS\COG_UNC\PLEIO_&trta.\&trtb._fold.tiff";


proc sgplot data=dsd.&trtb._fold_subs;where fold LE 10;
	XAXIS label = "Nominal -log10(pCOG)" min=0 max=3.2;
	YAXIS label = "Fold enrichment COG|CHOL";

	series x=lp_cog y=fold/group=c;
run;

ods printer printer=tiff file="W:\Final\RESULTS\COG_UNC\PLEIO_&trta.\&trtb._qq.tiff";

proc sgplot data=dsd.&trtb._fold_subs;where fold LE 10;
	XAXIS label = "Empirical -log10(qCOG|&trta.)" min=0 max=4;
	YAXIS label = "Nominal -log10(pCOG)" min=0 max=4;
	
	series x=lpqCOGsub y=lp_cog/group=c;
	series x=lpqCOGsub y=lpqCOGsub;
run;
ods printer close;
ods html;

	quit;

	%mend;

%mark(CHOL,chol);
%mark(LDL,ldl);
%mark(LOGCRP,logcrp);
%mark(LOGHDL,loghdl);


/** reverse  trait as primary and cognitive as secondary  ***/


%macro marka(trta,trtb);

libname dsd "W:\Final\RESULTS\COG_UNC\PLEIO_&trta.";
options nodate nonumber;
proc append base=dsd.&trtb._fold_allrev data=dsd.&trtb.rev2p0;run;
proc append base=dsd.&trtb._fold_allrev data=dsd.&trtb.rev2p5;run;
proc append base=dsd.&trtb._fold_allrev data=dsd.&trtb.rev3p0;run;
proc append base=dsd.&trtb._fold_allrev data=dsd.&trtb.rev3p5;run;
proc append base=dsd.&trtb._fold_allrev data=dsd.&trtb.rev4p0;run;


data dsd.&trtb._fold_subsrev;
	set dsd.&trtb._fold_allrev;
	c = cat("-log10(&trtb.addp) >= ",group);
	if fold > 10 then delete;
	/* if lp_TICS < 6.0 then output; */
run;

proc sort data=dsd.&trtb._fold_subsrev;
	by group lp_&trtb.;
run;


ods html close;
ods printer printer=tiff file="W:\Final\RESULTS\COG_UNC\PLEIO_&trta.\&trtb._foldrev.tiff";


proc sgplot data=dsd.&trtb._fold_subsrev;
	XAXIS label = "Nominal -log10(p&trta.)" min=0 max=4.1;
	YAXIS label = "Fold enrichment &trta.|COG";

	series x=lp_&trtb. y=fold/group=c;
run;

ods printer printer=tiff file="W:\Final\RESULTS\COG_UNC\PLEIO_&trta.\&trtb._qqrev.tiff";

proc sgplot data=dsd.&trtb._fold_subsrev;
	XAXIS label = "Empirical -log10(q&trta.|COG)" min=0 max=5;
	YAXIS label = "Nominal -log10(p&trta.)" min=0 max=5;
	
	series x=lpq&trtb.sub y=lp_&trtb./group=c;
	series x=lpq&trtb.sub y=lpq&trtb.sub;
run;
ods printer close;
ods html;

	quit;

	%mend;

%marka(CHOL,chol);
%marka(LDL,ldl);
%marka(LOGCRP,logcrp);
%marka(LOGHDL,loghdl);

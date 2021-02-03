
/*** Add FDR  **/

%macro mark(trta,trtb);

libname dsd "W:\Final\RESULTS\COG_UNC\PLEIO_&trta.";

data dsd.&trtb._fintab;
	set dsd.&trtb._fold_all;
	fdr = (10**-lp_cog) / (10**-lpqcogsub);
	*if lp_cog < 9 then output;
	run;

data ;set &trtb._fold_subs;
run;


data dsd.&trtb._fintabrev;
	set dsd.&trtb._fold_allrev;
	fdr = (10**-lp_&trtb.) / (10**-lpq&trtb.sub);
	*if lp_cog < 9 then output;
	run;

/***

proc sort data=dsd.crp_fold_subs;
	by SNP;
	run;
proc sort data=dsd.crp_impair_all;
	by SNP;
	run;
data dsd.crp_fintab;
	merge dsd.crp_fold_subs (in=insub) dsd.cog_impair_all (in=incog);
	by SNP;
	usedata = insub;
	run;
data dsd.crp_fintab;
	set dsd.crp_fintab;
	drop usedata;
	if usedata > 0 then output;
	run;
proc sort data=dsd.crp_fintab;
	by group lp_cog;
	run;

data mwl;
	set dsd.crp_fintab;
	if fold >= 2 then output;
	run;
proc print data=mwl;run;

proc sort data=dsd.crp_fintab;
	by fdr;
	run;
data mwl;
	set dsd.crp_fintab;
	if fdr <= 0.05 then output;
	run;
proc print data=mwl;
run;

****/

%mend;

%mark(CHOL,chol);
%mark(LDL ,ldl);
%mark(LOGCRP,logcrp);
%mark(LOGHDL,loghdl);


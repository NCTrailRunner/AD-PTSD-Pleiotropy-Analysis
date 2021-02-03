*libname dsd 'C:\Users\mwl17\Documents\2018 AD PTSD pleiotropy\AfAm GWAS';

libname dsd  "/isilon/datalake/irb/BG00316/Final/RESULTS/COG_UNC/PLEIO_CHOL";


/****
data dsd.cog_chol_logp;
	set dsd.cog_chol(rename=(addlogp=lp_chol));
	lp_cog = -log10(cogaddp); 
	*lp_PTSD = -log10(p_PTSD);
run;
data mwl;
	set dsd.cog_chol_logp;
proc sort data=mwl;
	by lp_cog;

data mwl;
	set mwl;
	rank = _N_;
run;
proc sort data=mwl;
	by lp_cog;

data mwl;
	set mwl nobs=ntot;
	by lp_cog;
	retain run_sum ctr;

	keep lp_cog run_sum ctr rank ntot ecdfy avg_rank lpqCOGall;
	
    if first.lp_cog=1 then do;
		run_sum = 0;
		ctr = 1;
		run_sum = rank;
		end;
	if first.lp_cog=0 and last.lp_cog=0 then do;
		run_sum = run_sum + rank;
		ctr = ctr + 1;
		end;
	if last.lp_cog=1 then do;
		if first.lp_cog = 0 then do;
			ctr = ctr + 1;
			run_sum = run_sum + rank;
			end;
		avg_rank = run_sum / ctr;
		*calculate the empirical CDF and lpq;
		ecdfy = avg_rank / (ntot + 1);
		lpqcogall = -log10(1 - ecdfy);
		output;
		end;

	proc sort data=mwl;
		by lp_cog;


	proc sort data=dsd.cog_chol_logp;
		by lp_cog;
		
		
	data dsd.pleio_all;
		merge dsd.cog_chol_logp mwl;
		by lp_cog;
run;

endsas;

***/


/*** Step 2: subsetting on secondary trait logp with threshold   ****/



%macro mark(thold,tholda);

data dsd.pleiosubs;
	set dsd.cog_chol_logp;
	*secondary phenotype inclusion threshold;
	if lp_chol >= &thold. then output;
run;
proc sort data=dsd.pleiosubs;
	by lp_cog;
*add rank variable;
data mwl;
	set dsd.pleiosubs;
	rank = _N_;
run;
proc sort data=mwl;
	by lp_cog;


*calculate average rank and empirical CDF;
data mwl;
	set mwl nobs=ntot;
	by lp_cog;
	retain run_sum ctr;
	*drop SNP lp_chol addp CHR pos;
	keep lp_cog run_sum ctr rank avg_rank ecdfy ntot lpqcogsub;
    if first.lp_cog=1 then do;
		run_sum = 0;
		ctr = 1;
		run_sum = rank;
		end;
	if first.lp_cog=0 & last.lp_cog=0 then do;
		run_sum = run_sum + rank;
		ctr = ctr + 1;
		end;
	if last.lp_cog=1 then do;
		if first.lp_cog = 0 then do;
			ctr = ctr + 1;
			run_sum = run_sum + rank;
			end;
		avg_rank = run_sum / ctr;
		ecdfy = avg_rank / (ntot + 1);
		lpqcogsub = -log10(1 - ecdfy);
		output;
		end;
*add in the new variables;
	proc sort data=dsd.pleiosubs;
		by lp_cog;
	proc sort data=mwl;
		by lp_cog;
	data dsd.pleio_sub;
		merge dsd.pleiosubs mwl;
		by lp_cog;
run;

proc sort data=dsd.pleio_sub;
	by SNP;
proc sort data=dsd.pleio_all;
	by SNP;
data dsd.pleiocomb;
	merge dsd.pleio_sub (in=insub) dsd.pleio_all (in=inall);
	by SNP;
	sub_idx = insub;
run;
data dsd.pleiocomb;
	set dsd.pleiocomb;
	if sub_idx = 1 then do;
	fold = 10 ** (lpqcogall - lpqcogsub);
	output;
	end;
run;


/***
proc gplot data=dsd.pleiocomb;
	plot fold*lp_AD/HAXIS= 0 to 12 by 1;
	plot (lpqADsub lp_AD)*lpqADsub/overlay;
run;

***/

data dsd.CHOL&tholda.;
	set dsd.pleiocomb;
	group = &thold.;
run;

%mend;

%mark(2.0,2p0);
%mark(2.5,2p5);
%mark(3.0,3p0);
%mark(3.5,3p5);
%mark(4.0,4p0);

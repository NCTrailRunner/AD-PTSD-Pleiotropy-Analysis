libname dsd 'C:\Users\mwl17\Documents\2018 AD PTSD pleiotropy\PTSD GWAS';
ODS graphics on;
/* initial dataset with p values for both phenotypes */
/* transform p values to -log10 */
data dsd.pleio_test_noAPOE;
	set dsd.pleio_test_noAPOE;
	lp_PTSD = -log10(p_PTSD); 
	lp_AD = -log10(p_AD);
run;
data mwl;
	set dsd.pleio_test_noAPOE;
proc sort data=mwl;
	by lp_PTSD;
/* add rank variable */
data mwl;
	set mwl;
	rank = _N_;
run;
proc sort data=mwl;
	by lp_PTSD;
/* calculate average rank and empirical CDF */
/* average rank is obtained by averaging the ranks of observations that have the same */
/* -log10(p PTSD).  Eventually, the average rank is normalized to N + 1 observations */
/* This algorithm matches the one used in JMP for save prob scores */
data mwl;
	set mwl nobs=ntot;
	by lp_PTSD;
	retain run_sum ctr;
	/* make sure to not overwrite these variables from original dataset */
	/* otherwise data gets scrambled in the merge operations */
	drop SNP lp_AD p_AD;
    if first.lp_PTSD=1 then do;
		run_sum = 0;
		ctr = 1;
		run_sum = rank;
		end;
	if first.lp_PTSD=0 & last.lp_PTSD=0 then do;
		run_sum = run_sum + rank;
		ctr = ctr + 1;
		end;
	if last.lp_PTSD=1 then do;
		if first.lp_PTSD = 0 then do;
			ctr = ctr + 1;
			run_sum = run_sum + rank;
			end;
		avg_rank = run_sum / ctr;
		/* calculate the empirical CDF and lpq */
		ecdfy = avg_rank / (ntot + 1);
		lpqPTSDall = - log10 (1 - ecdfy);
		output;
		end;
/* add in the new variables */
	proc sort data=dsd.pleio_test_noAPOE;
		by lp_PTSD;
	proc sort data=mwl;
		by lp_PTSD;
	data dsd.pleio_all;
		merge dsd.pleio_test_noAPOE mwl;
		by lp_PTSD;
run;
/* analysis of the subset */
data dsd.pleiosubs;
	set dsd.pleio_test_noAPOE;
	/* secondary phenotype inclusion threshold */
	if lp_AD >= 6.0 then output;
run;
proc sort data=dsd.pleiosubs;
	by lp_PTSD;
/* add rank variable */
data mwl;
	set dsd.pleiosubs;
	rank = _N_;
run;
proc sort data=mwl;
	by lp_PTSD;
/* calculate average rank and empirical CDF */
data mwl;
	set mwl nobs=ntot;
	by lp_PTSD;
	retain run_sum ctr;
	drop SNP lp_AD p_AD;
    if first.lp_PTSD=1 then do;
		run_sum = 0;
		ctr = 1;
		run_sum = rank;
		end;
	if first.lp_PTSD=0 & last.lp_PTSD=0 then do;
		run_sum = run_sum + rank;
		ctr = ctr + 1;
		end;
	if last.lp_PTSD=1 then do;
		if first.lp_PTSD = 0 then do;
			ctr = ctr + 1;
			run_sum = run_sum + rank;
			end;
		avg_rank = run_sum / ctr;
		ecdfy = avg_rank / (ntot + 1);
		lpqPTSDsub = - log10 (1 - ecdfy);
		output;
		end;
/* add in the new variables */
	proc sort data=dsd.pleiosubs;
		by lp_PTSD;
	proc sort data=mwl;
		by lp_PTSD;
	data dsd.pleio_sub;
		merge dsd.pleiosubs mwl;
		by lp_PTSD;
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
	fold = 10 ** (lpqPTSDall - lpqPTSDsub);
	output;
	end;
run;
proc gplot data=dsd.pleiocomb;
	plot fold*lp_PTSD;
	plot (lpqPTSDsub lp_PTSD)*lpqPTSDsub/overlay;
run;
data dsd.ADnE4_6;
	set dsd.pleiocomb;
	group = 6.0;
run;
quit;

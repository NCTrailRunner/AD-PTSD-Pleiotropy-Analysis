libname dsd 'C:\Users\mwl17\Documents\2018 AD PTSD pleiotropy\PTSD GWAS';
data dsd.ptsd_fold_subs;
	set dsd.ptsd_fold_all;
	fdr = (10**-lp_AD) / (10**-lpqADsub);
	/* if lp_AD < 9 then output; */

data mwl2;
	set dsd.ptsd_fold_subs;
	length nSNP $12;
	drop SNP;
	nSNP = SNP;
data dsd.ptsd_fold_subs;
	set mwl2;
	drop nSNP;
	SNP = nSNP;

proc sort data=dsd.ptsd_fold_subs;
	by SNP;
proc sort data=dsd.IGAP_GWAS;
	by SNP;
data dsd.ptsd_fintab;
	merge dsd.ptsd_fold_subs (in=insub) dsd.IGAP_GWAS (in=inIGAP);
	by SNP;
	usedata = insub;
data dsd.ptsd_fintab;
	set dsd.ptsd_fintab;
	drop usedata;
	if usedata > 0 then output;
proc sort data=dsd.ptsd_fintab;
	by group lp_AD;
data mwl;
	set dsd.ptsd_fintab;
	/* if fold >= 2 then output; */
/* proc print data=mwl; */
proc sort data=dsd.ptsd_fintab;
	by fdr;
data mwl;
	set dsd.ptsd_fintab;
	if fdr <= 0.05 then output;
proc print data=mwl;
quit;

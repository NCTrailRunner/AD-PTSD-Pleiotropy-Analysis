libname dsd 'C:\Users\mwl17\Documents\2018 AD PTSD pleiotropy\PTSD GWAS';
options nodate nonumber;
/* proc append base=dsd.ptsd_fold_all data=dsd.ptsd1;
proc append base=dsd.ptsd_fold_all data=dsd.ptsd1p5;
proc append base=dsd.ptsd_fold_all data=dsd.ptsd2;
proc append base=dsd.ptsd_fold_all data=dsd.ptsd3; */
/* proc append base=dsd.ptsd_fold_all data=dsd.ptsd1; */


data dsd.ptsd_fold_subs;
	set dsd.ptsd_fold_all;
	c = cat("-log10(Pptsd) >= ",group);
	/* if lp_AD < 6.0 then output; */

proc sort data=dsd.ptsd_fold_subs;
	by group lp_AD;

ods html close;
ods printer printer=tiff file='C:\Users\mwl17\Documents\2018 AD PTSD pleiotropy\ptsd gwas\ptsd fold.tiff';


proc sgplot data=dsd.ptsd_fold_subs;
	XAXIS label = "Nominal -log10(P_AD)" min=0 max=5;
	YAXIS label = "Fold enrichment AD|ptsd";

	series x=lp_AD y=fold/group=c;
run;
ods printer printer=tiff file='C:\Users\mwl17\Documents\2018 AD PTSD pleiotropy\ptsd gwas\ptsd_qq.tiff';

proc sgplot data=dsd.ptsd_fold_subs;
	XAXIS label = "Empirical -log10(qAD|ptsd)" min=0 max=4;
	YAXIS label = "Nominal -log10(P_AD)" min=0 max=4;
	
	series x=lpqADsub y=lp_AD/group=c;
	series x=lpqADsub y=lpqADsub;
run;
ods printer close;
ods html;

	quit;

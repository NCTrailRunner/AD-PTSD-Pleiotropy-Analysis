# AD-PTSD-Pleiotropy-Analysis

February 3, 2021

This repository contains the SAS code to complete the genetic pleiotropy analysis described in Lutz et al., Shared genetic etiology underlying late-onset Alzheimer's disease and posttraumatic stress syndrome” (1).  Details of the statistical analysis method are described in the paper.  

The SAS code is used for the steps necessary to calculate conditional probabilities for Alzheimer’s disease conditional with PTSD at different significance levels.  The code also produces the plots included in the manuscript.  Note that the Pleio_calc.sas programs run the primary conditional probability analysis for each direction, Q(AD|PTSD) and Q(PTSD|AD).  The programs that produce the plots and tables must be modified to reflect the specific dataset names used for each direction and therefore run twice, once for AD|PTSD and once for PTSD|AD.

For studies that are using multiple phenotypes, Mark Brown at Wake Forest University School of Medicine created a version of these programs that work with macros and has kindly shared this code.  

Details on using and modifying the code for specific applications are described in the file: Genetic Pleiotropy SAS code documentation.pdf.

1.	Lutz MW, Luo S, Williamson DE, Chiba-Falek O. Shared genetic etiology underlying late-onset Alzheimer's disease and posttraumatic stress syndrome. Alzheimers Dement. 2020 Sep;16(9):1280-1292. doi: 10.1002/alz.12128. Epub 2020 Jun 26. PMID: 32588970; PMCID: PMC7769164.

# Datasets in wheelP

A number of Data objects are provided. These correspond to the experimental data that was analyzed,
and underlying numeric values for several figures in the associated manuscript.
Raw sequence data can be found in the appropriate repositories ([GSE102248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102248), [PRJEB22060](https://www.ebi.ac.uk/ena/data/view/PRJEB22060)).

## *In vitro* assays

* **All.filtered** Data frame of individual OD600 measurements for 440 strains that were
assayed *in vitro* and that passed quality control. Assays were performed in 96-well plates
and measrumentes were taken over 72 hours. Data frame columns are as follows:
	* Row: row letter in the 96-well plate.
	* Column: column number in the 96-well plate.
	* OD600: optical density measurement at 600nm.
	* rep: replicate number (1 or 2).
	* hrs: time-point in hours when the measurement was taken.
	* Strain: strain ID.
	* condition: phosphate and exudate condition (plusP, minusP, plus2minusP or minus2plusP).
	
* **Features** Data frame of *in vitro* growth features for 440 strains. All of these values
are calculated from the growth data in **All.filtered**. The data frame includes the
following columns:
	* Strain: strain ID.
	* minus2plusP, munisP, plus2minusP, plusP: Area under the curve (AUC) under each condition
	* MP\_L3M, M\_L3M, P\_L3M, PM\_L3M: mean OD600 in the last three measurements in each condtion.
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_MAX, M\_MAX, P\_MAX, PM\_MAX: maximum optical density achieved in each condtion
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_HMT, M\_HMT, P\_HMT, PM\_HMT: mean time to reach half the maximum density in each condtion
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_MGS, M\_MGS, P\_MGS, PM\_MGS: maximum growth rate achieved in each condtion
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).

## Binary plant-bacterium assays

	Elongation
	
	Pi

	binP.all
	
	binP
	
## Synthetic community assays



	Map.colonization

	
	Strain.auc

	Tax.colonization



	dge.wheel

	wheelP.full

	wheelP.mapsplit

	wheelP.rna

## Underlying figure data

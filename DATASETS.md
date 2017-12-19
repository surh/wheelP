# Datasets in wheelP

A number of Data objects are provided. These correspond to the experimental data that was analyzed,
and underlying numeric values for several figures in the associated manuscript.
Raw sequence data can be found in the appropriate repositories ([GSE102248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102248), [PRJEB22060](https://www.ebi.ac.uk/ena/data/view/PRJEB22060)).

## *In vitro* assays

* **All.filtered** Data frame of individual OD600 measurements for 440 strains that were
assayed *in vitro* and that passed quality control. Assays were performed in 96-well plates
and measrumentes were taken over 72 hours. Data frame columns are as follows:
	* Row (row letter in the 96-well plate),
	* Column (column number in the 96-well plate), OD600 (optical
density measurement at 600nm), rep (replicate number), hrs (time-point in hours when
the measurement was taken), Strain (strain ID), condition (phosphate and exudate condition).
	
* **Features** Data frame of *in vitro* growth features 

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

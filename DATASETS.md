# Datasets in wheelP

A number of Data objects are provided. These correspond to the experimental data that was analyzed,
and underlying numeric values for several figures in the associated manuscript.
Raw sequence data can be found in the appropriate repositories ([GSE102248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102248), [PRJEB22060](https://www.ebi.ac.uk/ena/data/view/PRJEB22060)).

## *In vitro* assays

* **All.filtered** Data frame of individual OD600 measurements for 440 strains that were
assayed *in vitro* and that passed quality control. Assays were performed in 96-well plates
and measurements were taken over 72 hours. Data frame columns are as follows:
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
	* minus2plusP, munisP, plus2minusP, plusP: log-transformed area under the growth curve (AUC)
	under each condition. The raw values are in the **Strain.auc** data frame.
	* MP\_L3M, M\_L3M, P\_L3M, PM\_L3M: mean OD600 in the last three measurements in each condtion.
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_MAX, M\_MAX, P\_MAX, PM\_MAX: maximum optical density achieved in each condtion
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_HMT, M\_HMT, P\_HMT, PM\_HMT: mean time to reach half the maximum density in each condtion
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_MGS, M\_MGS, P\_MGS, PM\_MGS: maximum growth rate achieved in each condtion
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* PLF\_L3M, MLF\_L3M: log2 ratio of the mean OD600 in the last three measurements between each
	exudate condition and its control (PLF = log2(minus2plusP/plusP), MLF = log2(plus2minusP/minusP).
	* PLF\_MAX, MLF\_MAX: log2 ratio of the maximum optical density achieved between each
	exudate condition and its control (PLF = log2(minus2plusP/plusP), MLF = log2(plus2minusP/minusP).
	* PLF\_HMT, MLF\_HMT: log2 ratio of the mean time to reach half the maximum density between each
	exudate condition and its control (PLF = log2(minus2plusP/plusP), MLF = log2(plus2minusP/minusP).
	* PLF\_MGS, MLF\_MGS: log2 ratio of maximum growth rate achieved between each
	exudate condition and its control (PLF = log2(minus2plusP/plusP), MLF = log2(plus2minusP/minusP).

* **Strain.auc** Raw area under the growth curve (AUC) for each of 440 tesed strains. Row names
indicate strain ID. Columns of the data frame are as follows:
	* minus2plusP, munisP, plus2minusP, plusP: Area under the growth curve (AUC) for each condition.

## Binary plant-bacterium assays

* **binP.all** Individual plant shoot phosphate content measurements for plant-bacterium binary
assays involving 194 bacterial strains. Data frame columns are as follows:
	* sample: sample ID.
	* Experiment: biological replicate batch ID.
	* Replicate: biological replicate for the corresponding strain.
	* Date: assay date.
	* Experimenter: ID of individual who took the measurement.
	* StarP: starting phosphate and sucrose conditions used for germination.
	* EndP: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Treatment: Either bacterial strain ID or "No Bacteria" for controls.
	* Pi\_content: plant shoot pi content measurement in (mmol Pi) / (mg FW), where FW stands for
	shoot fresh weight.
	* Bacteria: Indicating whether bacteria was added (+Bacteria) or not (No Bacteria).

* **binP** Results of testing for the effect of individual bacterial strains on plant shoot phosphate
accumulation in plant-bacterium binary assays. All results are derived from data in **binP.all**.
Columns of the data frame are as follows:
	* Strain: strain ID.
	* minusP\_100uM, minusP\_30uM, plusP\_100uM, plusP\_30uM: log(fold-change) in plant shoot
	pi content caused by bacteria with respect to no bacteria controls in each condition.
	* minusP\_100uM.qval, minusP\_30uM.qval, plusP\_100uM.qval, plusP\_30uM.qval: Benjamini-Hochberg
	corrected p-value of the effect of each bacterial strains on plant shoot pi content
	in each condition.
	* Mean: mean log(fold-change) of shoot pi content caused by bacteria across all condtions.
	* Group: group to which a strain belongs (positive, negative, indifferent or none).
	
## Synthetic community assays

* **Elongation** Main root elongation measurements for individual plants treated with synthetic
communities. Measurements where made from pictures in imageJ. Columns of the data frame are as follows:
	* Picture: picture ID.
	* Treatment: string indicating the combination of phosphate conditions and bacterial tretment.
	* Elongation: main root elongation measurement in cm.
	* Experiment: biological replicate batch ID.
	* Plate: petri dish ID.
	* StarP: starting phosphate and sucrose conditions used for germination.
	* EndP: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Bacteria: ID of synthetic community added, or none.

* **Pi** Plate level phenptypic measurement for plants treated with synthetic communities. Columns
of the data frame are as follows:
	* id: measurement id
	* empty\_tubes: weight (mg) of empty eppendorf tubes prior to sample collection.
	* full\_tubes: weight (mg) of eppendorf tubes with sample.
	* mgFW: estimated sample weight (mg).
	* OD820nm: optical density measurement at 820nm used to estimate shoot Pi content.
	* Pi\_content: shoot phosphate content  in (mmol Pi) / (mg FW).
	* Experiment: biological replicate batch ID.
	* shoot\_area: Shoot area (cm^2) estimated with winrhizo.
	* Nplants: number of plants in plate.
	* Elongation: mean main root elongation for plants in the same plate (derived from
	**Elongation** data frame).
	* StarP: starting phosphate and sucrose conditions used for germination.
	* EndP: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Bacteria: ID of synthetic community added, or none.
	* order: sample processing order.
	* total.root: total root network length per plate (cm).
	* Plate: petri dish ID.

* **Map.colonization**: Metadata for samples that underwent microbial profiling via 16S gene
sequencing. The columns of the data frame are as follows:
	* ID: sample ID.
	* Well: position of the well that contained the sample for DNA extraction and all downstream
	sequencing library preparation.
	* SampleID: internal ID.
	* Genotype: plant genotype.
	* Bacteria: ID of synthetic community in sample or "No Bacteria".
	* Phosphate: descriptive string indicating phosphate treatment.
	* Replicate: replicate number within experimental batch.
	* Experiment: biological replicate batch for the current sample.
	* Plate: ID of plate for DNA extraction and all downstream library preparation steps.
	* Run: MiSeq run ID for 16S sequencing.
	* Barcode2: ID of the inner barcode used for multiplexing.
	* Frameshift: ID of the frameshifted primer combination used for library preparation.
	* Sample\_ID: ID for the sequencing machine.
	* Sample\_Name: sample name for the sequencing machine.
	* Sample\_plate: sample plate for the sequecing machine.
	* Sample\_Well: sample well in 96-well plate.
	* I7\_index\_ID: ID of index (outer barcode) used for multiplexing.
	* index: sequenc of index (outer barcode) used for multiplexing.
	* Fraction: sample fraction (Root, Agar or Inoculum).
	* rnaID: ID of corresponding RNA-seq library from the same sample.
	* AgarGroup: Grouping factor indicating samples that came from the same plate (agar environment).
	* P1, P2, P3, I1, I2, I3, N1, N2, N3: indicator variables showing whether a sample
	was treated with each of the 9 bacterial functional blocks.
	* Pre.Pi: starting phosphate and sucrose conditions used for germination.
	* Pos.Pi: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Inoculated: Inidicates whether bacteria was applied (+Bacteria) or not (No Bacteria).

* **Tax.colonization** Taxonomy and block allocation of bacteria used in synthetic community experiments.
The columns of the data frame are as follows:
	* ID: strain ID.
	* Taxonomy: taxonomy string.
	* Type: functional class of strain.
	* Block: functional block of strain.
	* MonoID: ID of strain in plant-bacterium binary assays.

* **dge.wheel** DGEList object. Full RNA-seq gene counts for synthetic community assays.
The count matrix is also available at GEO together with the raw sequence data 
([GSE102248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102248)). See the documentation
for the [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) package for
more details on the object.

* **wheelP.full** A Dataset object. Contains the sequence counts of each strain on each sample that
was profiled via 16S gene sequencing, as well as the sample metadata and strain taxonomic information.
The counts were obtained by mapping all reads against all reference sequences. See the documentation
of the [AMOR](https://github.com/surh/AMOR) package for more details on the object.

* **wheelP.mapsplit** A Dataset object. Contains the sequence counts of each strain on each sample that
was profiled via 16S gene sequencing, as well as the sample metadata and strain taxonomic information.
The counts were obtained by mapping reads against reference sequences of strains added to each specific
samples. See the documentation of the [AMOR](https://github.com/surh/AMOR) package for more
details on the object.

* **wheelP.rna** A Dataset object. Contains the sequence counts of each plant gene for each sample that
was profiled via RNA-seq, as well as the sample metadata and strain taxonomic information.
The counts were correspond to the counts in the **dge.wheel** object. See the documentation
of the [AMOR](https://github.com/surh/AMOR) package for more details on the object.

## Data from figures

Below we describe some data objects that contain the numeric values that underlie several
figures in the associated manuscript. All the numbers below were calculated from data in
the above objects and is therefore redundant but it is provided for the sake of completeness.

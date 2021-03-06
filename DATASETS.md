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
	* minus2plusP, minusP, plus2minusP, plusP: log-transformed area under the growth curve (AUC)
	under each condition. The raw values are in the **Strain.auc** data frame.
	* MP\_L3M, M\_L3M, P\_L3M, PM\_L3M: mean OD600 in the last three measurements in each condition.
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_MAX, M\_MAX, P\_MAX, PM\_MAX: maximum optical density achieved in each condition
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_HMT, M\_HMT, P\_HMT, PM\_HMT: mean time to reach half the maximum density in each condition
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* MP\_MGS, M\_MGS, P\_MGS, PM\_MGS: maximum growth rate achieved in each condition
	(MP = minus2plusP, M = minusP, P = plusP, PM = plus2minusP).
	* PLF\_L3M, MLF\_L3M: log2 ratio of the mean OD600 in the last three measurements between each
	exudate condition and its control (PLF = log2(minus2plusP/plusP), MLF = log2(plus2minusP/minusP).
	* PLF\_MAX, MLF\_MAX: log2 ratio of the maximum optical density achieved between each
	exudate condition and its control (PLF = log2(minus2plusP/plusP), MLF = log2(plus2minusP/minusP).
	* PLF\_HMT, MLF\_HMT: log2 ratio of the mean time to reach half the maximum density between each
	exudate condition and its control (PLF = log2(minus2plusP/plusP), MLF = log2(plus2minusP/minusP).
	* PLF\_MGS, MLF\_MGS: log2 ratio of maximum growth rate achieved between each
	exudate condition and its control (PLF = log2(minus2plusP/plusP), MLF = log2(plus2minusP/minusP).

* **Strain.auc** Raw area under the growth curve (AUC) for each of 440 tested strains. Row names
indicate strain ID. Columns of the data frame are as follows:
	* minus2plusP, minusP, plus2minusP, plusP: Area under the growth curve (AUC) for each condition.

* **gc.tree** A phylo object that defines a phylogenetic tree of 395 strains that were tested *in vitro*,
passed quality control, and ahad a full length 16S gene Sanger sequence.

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
	* Mean: mean log(fold-change) of shoot pi content caused by bacteria across all conditions.
	* Group: group to which a strain belongs (positive, negative, indifferent or none).
	
## Synthetic community assays

* **Elongation** Main root elongation measurements for individual plants treated with synthetic
communities. Measurements where made from pictures in imageJ. Columns of the data frame are as follows:
	* Picture: picture ID.
	* Treatment: string indicating the combination of phosphate conditions and bacterial treatment.
	* Elongation: main root elongation measurement in cm.
	* Experiment: biological replicate batch ID.
	* Plate: petri dish ID.
	* StarP: starting phosphate and sucrose conditions used for germination.
	* EndP: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Bacteria: ID of synthetic community added, or none.

* **Pi** Plate level phenotypic measurement for plants treated with synthetic communities. Columns
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
	* Sample\_plate: sample plate for the sequencing machine.
	* Sample\_Well: sample well in 96-well plate.
	* I7\_index\_ID: ID of index (outer barcode) used for multiplexing.
	* index: sequence of index (outer barcode) used for multiplexing.
	* Fraction: sample fraction (Root, Agar or Inoculum).
	* rnaID: ID of corresponding RNA-seq library from the same sample.
	* AgarGroup: Grouping factor indicating samples that came from the same plate (agar environment).
	* P1, P2, P3, I1, I2, I3, N1, N2, N3: indicator variables showing whether a sample
	was treated with each of the 9 bacterial functional blocks.
	* Pre.Pi: starting phosphate and sucrose conditions used for germination.
	* Pos.Pi: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Inoculated: indicates whether bacteria was applied (+Bacteria) or not (No Bacteria).

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

* **m1.wheel** A DGEGLM object. A fitted model of the RNA-seq data in **dge.wheel**.
See the documentation for the [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)
package for more details on the object.

## Data from figures

Below we describe some data objects that contain the numeric values that underlie several
figures in the associated manuscript. All the numbers below were calculated from data in
the above objects and is therefore redundant but it is provided for the sake of completeness.

* **cross.validation.error** These are the numeric values underlying figure 7B. The numbers
are the *leave-synthetic community-out* cross validation. For more details on its calculation
please visit the sister repository [wheelPi](https://github.com/clingsz/wheelPi).
The columns of the data frame are as follows:
	* Fold: the ID of the synthetic community removed
	* LM: The cross-validated error for the linear model in the corresponding fold.
	* INT: The cross-validated error for the linear model with interactions in the
	corresponding fold.
	* NN: The cross-validated error for the neural network in the corresponding fold.

* **group.colonization** These are the numeric values underlying figures S6B-C and fig12. They were
obtained directly from **wheelP.mapsplit** Dataset object. Columns of the data frame are as
follows:
	* Sample: internal sample ID.
	* ID: sequence library ID.
	* Bacteria: ID of syncom that was added to sample.
	* Fraction: Grouping factor indicating the sample fraction (Agar, Root or Inoculum), or
	whether the values correspond to theoretical expecttations ("Theoretical").
	* Pre.Pi: starting phosphate and sucrose conditions used for germination.
	* Pos.Pi: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Taxon: Grouping factor for abundances, corresponds to color bars in figures. It is either
	a taxonomic string at the Order level, or a functional block ID.
	* Abundance: Total count of strains of the given taxon or block in the specified sample.
	* Set: Indicates whether the sample is part of the original set of synthetic communities (Wheel),
	or part of the validation ("Validation") set.
	* Color: Inidicates whether the Taxon column corresponds to a taxonomic string for bacterial
	Order (Order) or functional block ID (Block).

* **mds.colonization** These are the numeric values that underlie figure S6A. They were calculated
directly from Dataset **wheelP.mapsplit**. Columns of the data frame are as follows.
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
	* Sample\_plate: sample plate for the sequencing machine.
	* Sample\_Well: sample well in 96-well plate.
	* I7\_index\_ID: ID of index (outer barcode) used for multiplexing.
	* index: sequence of index (outer barcode) used for multiplexing.
	* Fraction: sample fraction (Root, Agar or Inoculum).
	* rnaID: ID of corresponding RNA-seq library from the same sample.
	* AgarGroup: Grouping factor indicating samples that came from the same plate (agar environment).
	* P1, P2, P3, I1, I2, I3, N1, N2, N3: indicator variables showing whether a sample
	was treated with each of the 9 bacterial functional blocks.
	* Pre.Pi: starting phosphate and sucrose conditions used for germination.
	* Pos.Pi: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Inoculated: indicates whether bacteria was applied (+Bacteria) or not (No Bacteria).
	* Depth: number of reads in the samples.
	* Usable: number of non plant-derived reads in sample.
	* MDS1, MDS2, MDS3, MDS4, MDS5, MDS6: Multi-dimensional scaling axes.

* **nn.sensitivity** These are the numeric values underlying figure 7C. The numbers
are the results of the sensitivity analysis on the different models. For more details on its calculation
please visit the sister repository [wheelPi](https://github.com/clingsz/wheelPi).
The columns of the data frame are as follows:
	* ContextID: ID of the context (i.e. combination of input variables).
	* minusP.2.plusP: The effect of changing the starting phosphate contions from -Pi, to +Pi
	on plant shoot phosphate accumulation.
	* X30uM.2.100uM: The effect of changing the ending phosphate contions from 30 uM to 100 uM
	on plant shoot phosphate accumulation.
	* P1, P2, P3, I1, I2, I3, N1, N2, N3: The effect of each functional bacterial block on
	plant shoot phosphate accumulation.
	* Model: The model that generated the estimate (LM, INT or NN).

* **phen.additivity** Numeric values underlying figure 4 and figure 5B. Calculated directly from
the **Elongation** and **Pi** data frames. Columns of the data frame are as follows:
	* SynCom: Synthetic community ID.
	* StartP: starting phosphate and sucrose conditions used for germination.
	* EndP: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Estimate: linear model coefficient indicating the expected synthetic community effect
	on the corresponding plant phenotype.
	* SE: standard error of the estimate.
	* t.value: t-value of the estimate.
	* p.value: two-sided p.value testing the null hypothesis that the true value of
	the estimate is zero.
	* Measured: measured change in plant phenotype due to bacteria.
	* Predicted: predicted change in plant phenotype based on simple additive model.
	* SE.pred: standard error of the predicted phenotypic change.
	* Phenotype: phenotype being analyzed.

* **phen.presence.abundance** Numeric values underlying figure S7. Calculated directly from
the **Elongation** and **Pi** data frames, together with the **wheelP.mapsplit** Dataset.
Columns of the data frame are as follows:
	* SynCom: Synthetic community ID.
	* StartP: starting phosphate and sucrose conditions used for germination.
	* EndP: ending phosphate and sucrose conditions that were applied concomitant with bacteria.
	* Estimate: linear model coefficient indicating the expected synthetic community effect
	on the corresponding plant phenotype.
	* SE: standard error of the estimate.
	* t.value: t-value of the estimate.
	* p.value: two-sided p.value testing the null hypothesis that the true value of
	the estimate is zero.
	* Measured: measured change in plant phenotype due to bacteria.
	* Predicted: predicted change in plant phenotype based on simple additive model.
	* SE.pred: standard error of the predicted phenotypic change.
	* Type: Indicates whether the additive model included relative abundances (Abundance) or
	simply presence/absence (Block).
	* Phenotype: phenotype being analyzed.

* **pi.phu.cfu** Numeric values that underlie figure 2D and figure S3A. Directly obtained from
the experimenter. Columns of the data frame are as follows:
	* ID: sample ID.
	* Pre.treatment: Pre treatment utilized for germination: +Phi(1 mM Phi),
	+Pi(1 mM Pi), or -Pi(5 uM Pi).
	* Treatment: Pi level at the time of bacterial application.
	* Strain: Strain ID used in the sample.
	* Fraction: Sample fraction: SHOOT, AGAR, or ROOT.
	* CFU: colony forming units (c.f.u.) per fresh weight. Units are c.f.u./mL\*FW.
	* log.cfu: log of the c.f.u. per mL.
	* Replicate: biological replicate batch ID.

* **prediction.error** Theser are the numeric values underlying figure 7F. The numbers are
the mean different prediction error for the validation experiments and each model. 
For more details on its calculation please visit the sister repository
[wheelPi](https://github.com/clingsz/wheelPi). The columns of the data frame are as follows:
	* LM: The prediction error for the linear model.
	* INT: The prediction error for the linear model with interactions.
	* NN: The prediction error for the neural network.

* **pre.treatments** Numeric values underlying figure S3B. Directly obtained from the experimenter.
Columns of the data frame are as follows:
	* Sample: sample ID.
	* Experiment: biological replicate batch ID.
	* Pre.pi: Pre-treatment phosphate and phosphite levels used for germination.
	* Pi\_content: Plant shoot pi content. Units are mmol/mgFW.

* **signal.noise.ratio** These are the numeric values underlying figure S11. The numbers
are the signal and noise variances for each plant phenotype. The ratio of these variances 
(i.e. the signal to noise ratio) determines the feasibility of predictive modelling.
For more details on its calculation please visit the sister repository
[wheelPi](https://github.com/clingsz/wheelPi). The columns of the data frame are as follows:
	* Phenotype: The plant phenotype.
	* Which: Indicates whether the row corresponds to a noise (Noise.var) or signal (Signal.var)
	variance.
	* Variance: The signal or noise variance value for the corresponding phenotype.

* **validation.predicted.observed.** These are the numeric values underlying figure 7E. The numbers
are the predicted and observed plant shoot phosphate accumulation upon synthetic community
block replacements. For more details on its calculation please visit the sister repository
[wheelPi](https://github.com/clingsz/wheelPi). The columns of the data frame are as follows:
	* Change.ID: Descriptive string indicating the block replacement.
	* NN.pred: The predicted change in plant shoot phosphate accumulation by
	the neural network.
	* observed: The observed change in plant shoot phosphate accumulation.
	* Significant: flag indicating if the observed change in plant shoot phosphate
	accumulation was significant (1) or not (0).

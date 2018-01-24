# wheelP

Code and data for analysis described in the manuscript: "Design of synthetic bacterial
communities for predicable plant phenotypes".

This is an R package that is distributed with the aim on ensuring reproducibility.
The code was designed and tested to work specifically with the data in the package.

For the code that fits the neural network visit our sister repository
[wheelPi](https://github.com/clingsz/wheelPi).

# Installation

You will require the devtools package. Once that package is installed just type:

```r
devtools::install_github("surh/wheelP")
```

# Data

The raw sequence data is available in the appropriate academic repositories. RNA-seq data is at
NCBI GEO database under accession [GSE102248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102248).
Microbial 16S gene sequencing is available at the EBI SRA under accession [PRJEB22060](https://www.ebi.ac.uk/ena/data/view/PRJEB22060).

All other experimental data and numeric values underlying figures in the associated manuscript
is made available here as R package data. It can be accessed after installation using
`data(<DATASETNAME>)` in the R console. A full description of all the available
datasets can be found at the [DATASETS.md](DATASETS.md) file.

# Scripts

The [scripts directory](inst/scripts/) contains scripts used for analysis in the associated manuscript.
The scripts follow a naming convention in which a prefix is used to indicate the general set of
data that they analyze, and they are numbered to indicate they way in which they were run. In some
cases, scripts with higher numbers are dependent on output from scripts with lower numbers. The prefixes
are as follow:
* *gc\_* corresponds to the analysis of *in vitro* bacterial growth curves.
* *binP\_* corresponds to the analysis of plant phenotypes in plant-bacterium binary interaction assays.
* *syncom\_* corresponds to the analysis of plant phenotypes in bacterial synthetic community assays.
* *colonization\_* corresponds to the analysis of bacterial abundances in synthetic community assays.
* *rna\_* corresponds to the analysis of plant transcriptomes in synthetic community assays.

# Referencing

If using the code or data, please reference this repository's URL, as well as the following:

\*Herrera Paredes S, \*Gao T, Law TF, Finkel OM, Mucyn T, Texeira PJPL, Salas González I,
Feltcher ME, Powers MJ, Shank EA, Jones CD, Jojic V, Dangl JL & Castrillo G. "Design of 
synthetic bacterial communities for predicable plant phenotypes" (2017). *In press.
\*Co-first authors.

# Copyright & license

    (C) Copyright 2017 Sur Herrera Paredes

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


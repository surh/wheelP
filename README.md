# wheelP

Code and data for analysis of the effect of bacterial blocks on plant phosphate
accumulation and other related plant phenotypes.

This is an R package that is distributed with the aim on ensuring reproducibility.
The code was designed and tested to work specifically with the data in the package.

# Installation

You will require the devtools package. Once that package is installed just type:

```r
devtools::install_github("surh/wheelP")
```

# Data

The raw sequence data is available in the appropriate repositories. RNA-seq data is at NCBI GEO
database unde accession [GSE102248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102248).
Microbial profiling 16S gene sequencing is available at the EBI SRA under accession [PRJEB22060](https://www.ebi.ac.uk/ena/data/view/PRJEB22060).

Other experimental data and data underlying figures is made available here as R package data. It can be accessed after installation using `data(<DATASETNAME>)` in the R console. A full description of all the available datasets can be found at the [DATASETS.md](DATASETS.md) file.

# Referencing

If using the code of data, please reference this repository's URL, as well as the following:

Herrera Paredes S, Gao T, Law TF, Finkel OM, Mucyn T, Texeira PJPL, Salas González I,
Feltcher ME, Powers MJ, Shank EA, Jones CD, Jojic V, Dangl JL & Castrillo G. "A simplified
framework for dissecting complex host-microbiota interactions" (2017). *In revision*.

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


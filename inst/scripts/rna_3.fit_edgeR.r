# (C) Copyright 2017 Sur Herrera Paredes
# 
# This file is part of wheelP.
# 
# wheelP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# wheelP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with wheelP.  If not, see <http://www.gnu.org/licenses/>.

library(edgeR)
library(AMOR)
library(wheelP)
date()

# setwd("~/rhizogenomics/experiments/2017/today5//")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

# Load data
data(wheelP.rna)
Dat <- wheelP.rna
rm(wheelP.rna)

# Filter data for quick tests, comment out for real run.
# same for save statements which were used for quick test
# Dat <- measurable_taxa(Dat,min_reads_otu = 100,min_samples_otu = 400)

# Filter low abundant genes
Dat <- measurable_taxa(Dat,min_reads_otu = 5,min_samples_otu = 25)
Dat <- clean(Dat)

# Set formula
f1 <- formula(~ P1 + P2 + P3 +
                I1 + I2 + I3 +
                N1 + N2 + N3 +
                Phosphate +
                P1*Phosphate +
                P2*Phosphate +
                P3*Phosphate +
                I1*Phosphate +
                I2*Phosphate +
                I3*Phosphate +
                N1*Phosphate +
                N2*Phosphate +
                N3*Phosphate +
                Experiment)
# summary(lm(formula(paste(c("A260.280",paste(f1)), collapse = "")), Dat$Map))

# Set design matrix
design <- model.matrix(f1, data = Dat$Map)
dge <- DGEList(counts = Dat$Tab)

# Calculate normalization factors and dispersion
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge,design = design)
dge.wheel <- dge
# save(dge.wheel,file = "dge.test.rdat")
devtools::use_data(dge.wheel,
                   pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)
rm(dge.wheel)

# Fit Quasi-likelihood model
m1 <- glmQLFit(dge,design = design)
m1.wheel <- m1
# save(m1.wheel,file = "m1.test.rdat")
devtools::use_data(m1.wheel,
                   pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)
rm(m1.wheel)

date()

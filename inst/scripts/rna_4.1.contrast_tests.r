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

# Created this to add a test of positive vs negative at 30uM
library(edgeR)
library(wheelP)

date()

# Load previously fit model
data(m1.wheel)
m1 <- m1.wheel
rm(m1.wheel)

# Compare positive vs negative at 30uM
contrast <- c(0,1,1,1,0,0,0,-1,-1,-1,
              0,0,0,0,0,1,1,0,1,1,
              0,1,1,0,0,0,0,0,0,0,
              0,0,0,-1,-1,0,-1,-1,0,-1,
              -1)
qlf <- glmQLFTest(m1,contrast = contrast)
qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
qlf.tt <- qlf.tt$table
qlf.tt$Gene <- row.names(qlf.tt)
row.names(qlf.tt) <- NULL
hist(qlf.tt$PValue)
write.table(qlf.tt, "positive30uM_vs_negative30uM.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

date()

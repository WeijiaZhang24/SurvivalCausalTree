# SurvivalCausalTree

This repository implemented the survival causal tree (SCT) method proposed in Bioinformatics paper. The main computation part is written in C which significantly improves efficiency. The old code is at http://nugget.unisa.edu.au/Thus/SCT.zip

This code is based on rpart package from CRAN, and Susan Athey's causalTree package at https://github.com/susanathey/causalTree.

In addition, the propensity score in causalTree package has to be the same for all samples, in this package we fixed this problem so that it can be specified as a vector.

Use the script in the "test" folder for a demo.
                          

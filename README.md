# An Internal Validity Index Based on Density-involved Distance

## Requirement

The source code is written by Matlab r2016a. Versions lower than Matlab r2012a have not been tested.

## Simplest Demo

- RUN Ncut_test.m to determine the optimal partition from varied partitions (produced by Jianbo Shi's Normalized cuts). 

`CVDD.m` includes Algorithm 1: CVDD in our paper.
`Ncut_test.m` as an example includes Algorithm 2: CVDD-OP in our paper.

'OP_CA' (in `Ncut_test.m`) shows the comparison. [Purity, CVDD, CVNN, WB, Silhouette, CH, DB, Dunn, S_Dbw, I]

### Parameters in CVDD

No need to tune.

### Datasets used

File `Datasets_all30` includes 10 non-spherical clusters, 10 spherical clusters and 10 classification datasets (real datasets) used in the experiments of our paper.


## Issues, Questions, etc

Please report issues here on the github page or contact "hly4ml@gmail.com"

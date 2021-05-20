# An Internal Validity Index Based on Density-Involved Distance
https://ieeexplore.ieee.org/document/8672850

[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/an-internal-validity-index-based-on-density/clustering-algorithms-evaluation-on-pathbased)](https://paperswithcode.com/sota/clustering-algorithms-evaluation-on-pathbased?p=an-internal-validity-index-based-on-density)
[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/an-internal-validity-index-based-on-density/clustering-algorithms-evaluation-on-iris)](https://paperswithcode.com/sota/clustering-algorithms-evaluation-on-iris?p=an-internal-validity-index-based-on-density)
[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/an-internal-validity-index-based-on-density/clustering-algorithms-evaluation-on)](https://paperswithcode.com/sota/clustering-algorithms-evaluation-on?p=an-internal-validity-index-based-on-density)
[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/an-internal-validity-index-based-on-density/clustering-algorithms-evaluation-on-seeds)](https://paperswithcode.com/sota/clustering-algorithms-evaluation-on-seeds?p=an-internal-validity-index-based-on-density)
[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/an-internal-validity-index-based-on-density/clustering-algorithms-evaluation-on-97)](https://paperswithcode.com/sota/clustering-algorithms-evaluation-on-97?p=an-internal-validity-index-based-on-density)

## Requirement

The source code is written by Matlab r2016a. Versions lower than Matlab r2012a have not been tested.

## Simplest Demo

- RUN Ncut_test.m to determine the optimal partition from varied partitions (produced by Jianbo Shi's Normalized cuts). 

  `CVDD.m` includes Algorithm 1: CVDD in our paper.
  
  `Ncut_test.m` as an example includes Algorithm 2: CVDD-OP in our paper.

  'OP_CA' (in `Ncut_test.m`) shows the comparison. [ Purity, CVDD, CVNN, WB, Silhouette, CH, DB, Dunn, S_Dbw, I ]

### Parameters in CVDD

No need to tune.

### Datasets used

File `Datasets_all30` includes 10 non-spherical clusters, 10 spherical clusters and 10 classification datasets (real datasets) used in the experiments of our paper.


## Issues, Questions, etc

Please report issues here on the github page or contact "hly4ml@gmail.com"

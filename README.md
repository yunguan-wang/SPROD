# SPROD
SPROD: De-noising of Spatial Expression Profiling Data Based on Latent Graph Learning of in situ Position and Image Data

## Introduction
Spatial expression profiling (SEP) techniques provide gene expression close to or even superior to single cell resolution, while retaining the physical locations of sequencing and sometimes also provide matched pathological images. However, the expression data captured by spatial sequencing techniques suffer from severe inaccuracies, including but not limited to drop-outs as in regular single cell RNA-sequencing (scRNA-seq) data. To reduce the level of noise in the data, we developped the SPROD tool, which incorporated image information and spot positional information and used latent graph modeling to correct the gene expression data of each spot based on its neighboring spots.

## Graphical abstract

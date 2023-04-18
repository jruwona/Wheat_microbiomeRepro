---
title: "WheatMicrobiomeRepro_analayis"
author: "Justice Ruwona"
date: "2023-04-10"
output: 
  html_document:
      keep_md: yes
---
### Load in Libraries

```r
# DownLoad packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
```

```
## Warning: package 'BiocManager' was built under R version 4.2.3
```

```r
BiocManager::install("phyloseq")
```

```
## Bioconductor version 3.16 (BiocManager 1.30.20), R 4.2.2 (2022-10-31 ucrt)
```

```
## Warning: package(s) not installed when version(s) same as or greater than current; use
##   `force = TRUE` to re-install: 'phyloseq'
```

```
## Installation paths not writeable, unable to update packages
##   path: C:/Program Files/R/R-4.2.2/library
##   packages:
##     boot, class, codetools, foreign, lattice, MASS, Matrix, mgcv, nlme,
##     spatial, survival
```

```
## Old packages: 'cli', 'fastmap', 'gargle', 'htmltools', 'rhdf5', 'xfun', 'zoo'
```

```r
#load libraries
library(phyloseq)
library(tidyverse)
```

```
## Warning: package 'ggplot2' was built under R version 4.2.3
```

```
## Warning: package 'tibble' was built under R version 4.2.3
```

```
## Warning: package 'dplyr' was built under R version 4.2.3
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.1     ✔ readr     2.1.4
## ✔ forcats   1.0.0     ✔ stringr   1.5.0
## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
## ✔ purrr     1.0.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the ]8;;http://conflicted.r-lib.org/conflicted package]8;; to force all conflicts to become errors
```

```r
library(Biostrings)
```

```
## Loading required package: BiocGenerics
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:lubridate':
## 
##     intersect, setdiff, union
## 
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
## 
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which.max, which.min
## 
## Loading required package: S4Vectors
## Loading required package: stats4
## 
## Attaching package: 'S4Vectors'
## 
## The following objects are masked from 'package:lubridate':
## 
##     second, second<-
## 
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
## 
## The following object is masked from 'package:tidyr':
## 
##     expand
## 
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
## 
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## 
## The following object is masked from 'package:lubridate':
## 
##     %within%
## 
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
## 
## The following object is masked from 'package:purrr':
## 
##     reduce
## 
## The following object is masked from 'package:phyloseq':
## 
##     distance
## 
## The following object is masked from 'package:grDevices':
## 
##     windows
## 
## Loading required package: XVector
## 
## Attaching package: 'XVector'
## 
## The following object is masked from 'package:purrr':
## 
##     compact
## 
## Loading required package: GenomeInfoDb
## 
## Attaching package: 'Biostrings'
## 
## The following object is masked from 'package:base':
## 
##     strsplit
```

```r
library(ggplot2)
library(vegan)
```

```
## Warning: package 'vegan' was built under R version 4.2.3
```

```
## Loading required package: permute
```

```
## Warning: package 'permute' was built under R version 4.2.3
```

```
## Loading required package: lattice
```

```
## Warning: package 'lattice' was built under R version 4.2.3
```

```
## This is vegan 2.6-4
```

```r
library(dplyr)
library(lme4)
```

```
## Warning: package 'lme4' was built under R version 4.2.3
```

```
## Loading required package: Matrix
```

```
## Warning: package 'Matrix' was built under R version 4.2.3
```

```
## 
## Attaching package: 'Matrix'
## 
## The following object is masked from 'package:S4Vectors':
## 
##     expand
## 
## The following objects are masked from 'package:tidyr':
## 
##     expand, pack, unpack
```

```r
library(emmeans)
```

```
## Warning: package 'emmeans' was built under R version 4.2.3
```

```r
library(multcomp)
```

```
## Warning: package 'multcomp' was built under R version 4.2.3
```

```
## Loading required package: mvtnorm
## Loading required package: survival
```

```
## Warning: package 'survival' was built under R version 4.2.3
```

```
## Loading required package: TH.data
```

```
## Warning: package 'TH.data' was built under R version 4.2.3
```

```
## Loading required package: MASS
```

```
## Warning: package 'MASS' was built under R version 4.2.3
```

```
## 
## Attaching package: 'MASS'
## 
## The following object is masked from 'package:dplyr':
## 
##     select
## 
## 
## Attaching package: 'TH.data'
## 
## The following object is masked from 'package:MASS':
## 
##     geyser
```

```r
##### Set global options #####

# no scientific notation
options(scipen=10000) 

# color blind pallet used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
```
### READ IN DATA STEPS

```r
###### Read in data ####

#Loading the mapping file metadata
samp_dat <- read.csv("Data_files/mapping_its.csv")
samp_dat <- samp_dat %>%
  filter(Crop == "Yr_1") # filter only wheat dataset, exclude other crops in yr2, yr3 of rotation
rownames(samp_dat) <- samp_dat$Sample_Name # format row names must match OTU table headers
SAMP.fungi <- phyloseq::sample_data(samp_dat)

# OTU table - 
otu <- read.csv("Data_files/otu_table_its.csv")
treatment <- factor(c("T1", "T2", "T3", "T4"))
rownames(otu) <- otu$OTU_ID # format row names must match OTU table headers
otu <- otu[,-1] # shifting to 1, numeric matrix
OTU.fungi <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)

# Taxonomy
taxonomy.fungi <- read.csv("Data_files/taxonomy_table_its.csv")
rownames(taxonomy.fungi) <- taxonomy.fungi$OTU_ID
taxonomy.fungi2 <- taxonomy.fungi %>%
  subset(Kingdom == "Fungi")

TAX.fungi <- phyloseq::tax_table(as.matrix(taxonomy.fungi2))

# Fasta
#FASTA.fungi <- readDNAStringSet("ITS_all_crops.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
```

```r
###### Create Initial Phyloseq object ######
fungi.unedited <- phyloseq::phyloseq(OTU.fungi, TAX.fungi, SAMP.fungi)
str(fungi.unedited)
```

```
## Formal class 'phyloseq' [package "phyloseq"] with 5 slots
##   ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
##   .. .. ..@ .Data        : num [1:7233, 1:216] 437 1293 1035 15 2042 ...
##   .. .. .. ..- attr(*, "dimnames")=List of 2
##   .. .. .. .. ..$ : chr [1:7233] "F_OTU_20" "F_OTU_47" "F_OTU_1" "F_OTU_16" ...
##   .. .. .. .. ..$ : chr [1:216] "W_C1_T1_R1_L" "W_C1_T1_R1_R" "W_C1_T1_R1_S" "W_C1_T1_R2_L" ...
##   .. .. ..@ taxa_are_rows: logi TRUE
##   .. .. ..$ dim     : int [1:2] 7233 216
##   .. .. ..$ dimnames:List of 2
##   .. .. .. ..$ : chr [1:7233] "F_OTU_20" "F_OTU_47" "F_OTU_1" "F_OTU_16" ...
##   .. .. .. ..$ : chr [1:216] "W_C1_T1_R1_L" "W_C1_T1_R1_R" "W_C1_T1_R1_S" "W_C1_T1_R2_L" ...
##   ..@ tax_table:Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot
##   .. .. ..@ .Data: chr [1:7233, 1:8] "F_OTU_20" "F_OTU_47" "F_OTU_1" "F_OTU_16" ...
##   .. .. .. ..- attr(*, "dimnames")=List of 2
##   .. .. .. .. ..$ : chr [1:7233] "F_OTU_20" "F_OTU_47" "F_OTU_1" "F_OTU_16" ...
##   .. .. .. .. ..$ : chr [1:8] "OTU_ID" "Kingdom" "Phylum" "Class" ...
##   .. .. ..$ dim     : int [1:2] 7233 8
##   .. .. ..$ dimnames:List of 2
##   .. .. .. ..$ : chr [1:7233] "F_OTU_20" "F_OTU_47" "F_OTU_1" "F_OTU_16" ...
##   .. .. .. ..$ : chr [1:8] "OTU_ID" "Kingdom" "Phylum" "Class" ...
##   ..@ sam_data :'data.frame':	216 obs. of  12 variables:
## Formal class 'sample_data' [package "phyloseq"] with 4 slots
##   .. .. ..@ .Data    :List of 12
##   .. .. .. ..$ : chr [1:216] "W_C1_T1_R1_L" "W_C1_T1_R1_R" "W_C1_T1_R1_S" "W_C1_T1_R2_L" ...
##   .. .. .. ..$ : chr [1:216] "Yr_1" "Yr_1" "Yr_1" "Yr_1" ...
##   .. .. .. ..$ : chr [1:216] "C1" "C1" "C1" "C1" ...
##   .. .. .. ..$ : chr [1:216] "T1" "T1" "T1" "T1" ...
##   .. .. .. ..$ : chr [1:216] "R1" "R1" "R1" "R2" ...
##   .. .. .. ..$ : chr [1:216] "Leaf" "Root" "Stem" "Leaf" ...
##   .. .. .. ..$ : chr [1:216] "Wheat_C1" "Wheat_C1" "Wheat_C1" "Wheat_C1" ...
##   .. .. .. ..$ : chr [1:216] "Wheat_T1" "Wheat_T1" "Wheat_T1" "Wheat_T1" ...
##   .. .. .. ..$ : chr [1:216] "T1_Leaf" "T1_Root" "T1_Stem" "T1_Leaf" ...
##   .. .. .. ..$ : chr [1:216] "aboveground" "belowground" "aboveground" "aboveground" ...
##   .. .. .. ..$ : chr [1:216] "Kalamazoo" "Kalamazoo" "Kalamazoo" "Kalamazoo" ...
##   .. .. .. ..$ : chr [1:216] "chemical" "chemical" "chemical" "chemical" ...
##   .. .. ..@ names    : chr [1:12] "Sample_Name" "Crop" "Collection" "Treatment" ...
##   .. .. ..@ row.names: chr [1:216] "W_C1_T1_R1_L" "W_C1_T1_R1_R" "W_C1_T1_R1_S" "W_C1_T1_R2_L" ...
##   .. .. ..@ .S3Class : chr "data.frame"
##   ..@ phy_tree : NULL
##   ..@ refseq   : NULL
```




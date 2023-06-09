---
title: "WheatMicrobiomeRepro_analayis"
author: "Justice Ruwona"
date: "2023-04-10"
output: 
  html_document:
      keep_md: yes
---
# 1.0 Load in Libraries

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
library(ggpubr)
```

```
## Warning: package 'ggpubr' was built under R version 4.2.3
```

```r
##### Set global options #####

# no scientific notation
options(scipen=10000) 

# color blind pallet  
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
```
# 2.0 READ IN DATA STEPS

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

# Taxonomy table and choose Fungi only dataset
taxonomy.fungi <- read.csv("Data_files/taxonomy_table_its.csv")
rownames(taxonomy.fungi) <- taxonomy.fungi$OTU_ID
taxonomy.fungi2 <- taxonomy.fungi %>%
  subset(Kingdom == "Fungi")

TAX.fungi <- phyloseq::tax_table(as.matrix(taxonomy.fungi2))

# Fasta
#FASTA.fungi <- readDNAStringSet("ITS_all_crops.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
```
### Create Phyloseq object 

```r
###### Create Initial Phyloseq object ######
fungi.unedited <- phyloseq::phyloseq(OTU.fungi, TAX.fungi, SAMP.fungi)
```

### Prune sample reads and run reads statistics

```r
## Remove samples with less than 5000 reads
fungi.edited <- fungi.unedited %>% 
  prune_samples(sample_sums(.) > 5000, .) %>% # remove samples below 5,000 reads
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove taxa with less than 1 reads
sample_sums(fungi.edited) %>%
  sort()
```

```
## W_C1_T3_R2_R W_C1_T4_R6_R W_C1_T1_R4_L W_C2_T3_R2_R W_C2_T1_R3_S W_C1_T3_R2_L 
##         5180         5702         6805         6971         7374         7518 
## W_C3_T3_R2_R W_C2_T1_R3_R W_C1_T3_R3_R W_C3_T3_R5_S W_C1_T4_R5_R W_C1_T3_R4_R 
##         7858         7913         7996         8019         8298         8424 
## W_C1_T1_R2_S W_C2_T1_R2_S W_C2_T4_R4_R W_C1_T4_R3_L W_C1_T2_R6_S W_C3_T4_R4_R 
##         8918         9304         9367         9389         9451         9550 
## W_C1_T1_R3_R W_C2_T1_R5_L W_C3_T1_R6_R W_C3_T4_R1_R W_C1_T1_R1_L W_C2_T2_R6_L 
##         9731         9805        10236        10384        10454        10562 
## W_C1_T1_R5_L W_C2_T3_R5_S W_C1_T3_R1_L W_C1_T1_R1_R W_C1_T1_R6_R W_C1_T3_R6_S 
##        10584        10600        10697        10716        10779        10820 
## W_C2_T1_R2_R W_C3_T4_R2_L W_C1_T3_R3_S W_C3_T2_R3_S W_C1_T1_R1_S W_C3_T1_R1_S 
##        10963        11161        11201        11230        11343        11365 
## W_C2_T2_R3_S W_C3_T2_R1_R W_C1_T2_R4_S W_C1_T2_R5_R W_C1_T1_R6_S W_C1_T3_R5_S 
##        11383        11627        11794        11829        11919        11954 
## W_C3_T3_R4_R W_C2_T1_R1_R W_C2_T1_R6_R W_C3_T2_R6_L W_C1_T4_R6_S W_C1_T4_R4_R 
##        11957        12233        12247        12318        12531        12594 
## W_C2_T4_R2_L W_C1_T1_R6_L W_C1_T2_R3_L W_C1_T3_R5_R W_C3_T2_R6_S W_C2_T3_R5_R 
##        12598        12617        12727        12791        12819        12933 
## W_C3_T2_R5_S W_C2_T4_R5_R W_C1_T2_R1_L W_C1_T2_R1_R W_C1_T3_R5_L W_C2_T3_R3_L 
##        12950        12981        12993        13003        13036        13150 
## W_C2_T3_R1_L W_C1_T2_R3_S W_C3_T1_R1_R W_C3_T2_R4_S W_C1_T2_R3_R W_C1_T4_R2_S 
##        13387        13408        13571        13708        13782        13999 
## W_C3_T3_R3_L W_C2_T2_R1_L W_C2_T1_R4_L W_C1_T3_R3_L W_C1_T3_R4_S W_C3_T1_R2_S 
##        14117        14167        14220        14279        14339        14566 
## W_C2_T2_R3_L W_C1_T3_R2_S W_C2_T2_R3_R W_C2_T2_R6_S W_C3_T3_R3_R W_C1_T1_R4_S 
##        14577        14606        14695        14791        14823        14896 
## W_C3_T2_R3_R W_C1_T4_R1_R W_C1_T4_R4_L W_C1_T1_R3_S W_C3_T3_R2_S W_C2_T1_R6_L 
##        15011        15303        15360        15386        15424        15540 
## W_C1_T4_R4_S W_C2_T4_R5_L W_C1_T4_R2_L W_C1_T4_R3_R W_C3_T1_R6_S W_C2_T1_R6_S 
##        15575        15604        15607        15619        15720        15782 
## W_C2_T3_R5_L W_C1_T3_R4_L W_C2_T4_R4_S W_C3_T1_R4_L W_C1_T4_R5_L W_C3_T2_R1_L 
##        15845        15916        15986        16070        16088        16096 
## W_C1_T1_R2_L W_C3_T2_R3_L W_C1_T4_R3_S W_C2_T3_R3_R W_C2_T2_R4_S W_C2_T2_R5_R 
##        16180        16289        16337        16346        16390        16390 
## W_C3_T4_R4_S W_C3_T3_R1_L W_C2_T4_R3_L W_C2_T3_R2_S W_C2_T1_R1_S W_C1_T1_R5_S 
##        16420        16453        16512        16536        16722        16806 
## W_C2_T1_R5_R W_C1_T1_R4_R W_C2_T1_R4_S W_C2_T2_R2_L W_C2_T1_R4_R W_C1_T4_R6_L 
##        16813        17039        17167        17213        17343        17385 
## W_C3_T1_R4_S W_C1_T3_R6_R W_C3_T3_R5_L W_C1_T1_R5_R W_C3_T1_R3_R W_C1_T4_R1_S 
##        17416        17585        17617        18124        18188        18268 
## W_C3_T2_R2_R W_C1_T2_R4_L W_C1_T3_R1_S W_C1_T4_R2_R W_C1_T2_R2_L W_C2_T4_R6_R 
##        18443        18457        18602        18726        18798        18865 
## W_C2_T1_R2_L W_C2_T1_R1_L W_C2_T2_R1_S W_C1_T3_R1_R W_C1_T1_R2_R W_C1_T2_R2_R 
##        19074        19167        19386        19682        19758        19772 
## W_C2_T3_R6_S W_C2_T4_R4_L W_C1_T2_R4_R W_C2_T2_R6_R W_C3_T3_R5_R W_C2_T1_R5_S 
##        19815        19856        19865        20240        20576        20602 
## W_C3_T1_R1_L W_C1_T2_R5_S W_C2_T2_R5_S W_C1_T4_R5_S W_C2_T1_R3_L W_C1_T2_R1_S 
##        21224        21239        21254        21481        21553        21597 
## W_C1_T1_R3_L W_C3_T3_R3_S W_C2_T3_R4_R W_C3_T4_R3_L W_C1_T3_R6_L W_C1_T2_R5_L 
##        21929        22433        22764        23150        23307        23436 
## W_C3_T4_R2_S W_C2_T4_R2_R W_C3_T4_R1_S W_C3_T1_R5_L W_C2_T3_R6_L W_C2_T2_R4_L 
##        23838        25337        25842        26017        26186        26230 
## W_C3_T4_R6_R W_C3_T4_R4_L W_C2_T4_R5_S W_C1_T2_R2_S W_C3_T1_R3_S W_C3_T3_R4_L 
##        26804        26832        27550        27678        28109        28230 
## W_C2_T3_R2_L W_C3_T4_R5_L W_C3_T3_R6_L W_C3_T2_R1_S W_C2_T3_R4_S W_C3_T2_R2_S 
##        28283        28384        28507        28584        28970        29610 
## W_C3_T2_R6_R W_C1_T4_R1_L W_C3_T2_R4_L W_C3_T3_R4_S W_C2_T3_R3_S W_C2_T2_R2_R 
##        29743        29911        30042        30255        30510        30644 
## W_C3_T3_R2_L W_C3_T4_R5_R W_C3_T3_R1_R W_C2_T4_R1_S W_C3_T1_R2_R W_C3_T3_R6_R 
##        30831        31244        32308        32621        33360        33601 
## W_C3_T4_R5_S W_C3_T1_R2_L W_C2_T3_R1_R W_C2_T3_R1_S W_C2_T3_R4_L W_C2_T4_R3_S 
##        33824        34142        34681        36045        36465        36550 
## W_C2_T4_R3_R W_C2_T4_R6_L W_C3_T1_R4_R W_C3_T1_R6_L W_C3_T4_R3_S W_C3_T1_R3_L 
##        36563        37641        37808        38104        38570        38666 
## W_C2_T2_R4_R W_C3_T3_R1_S W_C3_T2_R5_L W_C2_T3_R6_R W_C3_T2_R2_L W_C3_T4_R3_R 
##        38724        39209        40904        41286        42124        42163 
## W_C3_T3_R6_S W_C2_T2_R2_S W_C2_T4_R6_S W_C3_T2_R4_R W_C3_T4_R6_L W_C2_T2_R5_L 
##        42256        42414        42932        43633        45053        45686 
## W_C2_T4_R1_L W_C3_T1_R5_R W_C3_T1_R5_S W_C3_T4_R1_L W_C3_T4_R6_S 
##        48152        49483        49950        53917        63519
```

```r
# New number of total reads
sum(sample_sums(fungi.edited)) # should read 4,272,434
```

```
## [1] 4272630
```

```r
# Mean and median read depth 
mean(sample_sums(fungi.edited)) #  should read 20,443.21
```

```
## [1] 20443.21
```

```r
median(sample_sums(fungi.edited)) #  should read 16,512
```

```
## [1] 16512
```

```r
# Histogram including median read depth
read.depths <- data.frame(sample_sums(fungi.edited))
colnames(read.depths) <- "read.depth"
read.depth.plot <- ggplot(read.depths, aes(read.depth)) +
  geom_histogram(fill = cbbPalette[[5]], color = "black") + 
  geom_vline(xintercept = median(sample_sums(fungi.unedited)), linetype = "dashed") + 
  theme_classic() + 
  ggtitle("Reads Depth Histogram") +
  xlab("Read Depth") + ylab("Count")
read.depth.plot
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Wheat_microbiomeReproScripts_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

# 3.0 DATA ANALYSIS

```r
# Global fungal community  dataset combining growth stage, plant organ, and crop management strategy 
subset_phylo <- subset_samples(fungi.edited)
```

## 3.1 Alpha diversity for global fungal endophytic communities

```r
# Rename fungi.edited data
physeq_sub <- fungi.unedited 
      
# Calculate alpha diversity metrics (Shannon, Simpson, and Observed OTUs) and print table
alpha_div <- estimate_richness(physeq_sub, measures=c("Shannon", "Simpson", "Observed"))
# Add treatment, collection time interval and tissue information to alpha diversity table
alpha_div <- cbind(Treatment = as.factor(sample_data(physeq_sub)$Treatment), 
                   Collection = as.factor(sample_data(physeq_sub)$Collection),
                   Tissue = as.factor(sample_data(physeq_sub)$Tissue),
                   alpha_div)
# Calculate summary statistics for each alpha diversity metric 
alpha_div_summary <- alpha_div %>%
  group_by(Treatment, Collection, Tissue) %>%
  summarize(mean_shannon = mean(Shannon), 
            sd_shannon = sd(Shannon), 
            mean_simpson = mean(Simpson), 
            sd_simpson = sd(Simpson), 
            mean_observed = mean(Observed), 
            sd_observed = sd(Observed))
```

```
## `summarise()` has grouped output by 'Treatment', 'Collection'. You can override
## using the `.groups` argument.
```

```r
alpha_div_summary # view summary of alpha diversity
```

```
## # A tibble: 36 × 9
## # Groups:   Treatment, Collection [12]
##    Treatment Collection Tissue mean_shannon sd_shannon mean_simpson sd_simpson
##    <fct>     <fct>      <fct>         <dbl>      <dbl>        <dbl>      <dbl>
##  1 T1        C1         Leaf           2.96      0.414        0.866     0.0670
##  2 T1        C1         Root           2.71      0.497        0.797     0.119 
##  3 T1        C1         Stem           2.46      0.360        0.783     0.105 
##  4 T1        C2         Leaf           2.87      0.376        0.869     0.0597
##  5 T1        C2         Root           3.30      0.420        0.923     0.0281
##  6 T1        C2         Stem           3.05      0.209        0.898     0.0237
##  7 T1        C3         Leaf           1.94      0.588        0.642     0.154 
##  8 T1        C3         Root           2.95      0.451        0.876     0.0580
##  9 T1        C3         Stem           2.34      0.631        0.698     0.186 
## 10 T2        C1         Leaf           2.69      1.37         0.920     0.0526
## # ℹ 26 more rows
## # ℹ 2 more variables: mean_observed <dbl>, sd_observed <dbl>
```

### Using whole plant tissue (Global community) Shannon index to compare if management communities are different

```r
#### ANOVA based on Shannon index
lm.alphadiv <- (lm(mean_shannon ~ Treatment*Collection,data = alpha_div_summary))
anova(lm.alphadiv)
```

```
## Analysis of Variance Table
## 
## Response: mean_shannon
##                      Df Sum Sq Mean Sq F value   Pr(>F)   
## Treatment             3 0.3975 0.13251  0.9673 0.424303   
## Collection            2 2.0945 1.04726  7.6450 0.002699 **
## Treatment:Collection  6 1.3416 0.22361  1.6323 0.181554   
## Residuals            24 3.2877 0.13699                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
### Test ls means using Tukey
smeans <- emmeans(lm.alphadiv, ~ Treatment | Collection) # estimate lsmeans of Treatment within Collection
Results_lsmeans <- cld(smeans, alpha = 0.05, reversed = TRUE, details = TRUE) # contrast with Tukey adjustment
Results_lsmeans # expect to see Collection 3 T2 - T4   p = 0.0233 is Significant
```

```
## $emmeans
## Collection = C1:
##  Treatment emmean    SE df lower.CL upper.CL .group
##  T4          2.97 0.214 24     2.53     3.42  1    
##  T3          2.76 0.214 24     2.32     3.20  1    
##  T1          2.71 0.214 24     2.27     3.15  1    
##  T2          2.65 0.214 24     2.21     3.09  1    
## 
## Collection = C2:
##  Treatment emmean    SE df lower.CL upper.CL .group
##  T3          3.23 0.214 24     2.79     3.67  1    
##  T2          3.23 0.214 24     2.78     3.67  1    
##  T4          3.19 0.214 24     2.75     3.63  1    
##  T1          3.08 0.214 24     2.64     3.52  1    
## 
## Collection = C3:
##  Treatment emmean    SE df lower.CL upper.CL .group
##  T2          3.10 0.214 24     2.66     3.54  1    
##  T3          2.76 0.214 24     2.32     3.20  12   
##  T1          2.41 0.214 24     1.97     2.85  12   
##  T4          2.16 0.214 24     1.72     2.60   2   
## 
## Confidence level used: 0.95 
## P value adjustment: tukey method for comparing a family of 4 estimates 
## significance level used: alpha = 0.05 
## NOTE: If two or more means share the same grouping symbol,
##       then we cannot show them to be different.
##       But we also did not show them to be the same. 
## 
## $comparisons
## Collection = C1:
##  contrast estimate    SE df t.ratio p.value
##  T1 - T2   0.05519 0.302 24   0.183  0.9978
##  T3 - T2   0.10381 0.302 24   0.344  0.9857
##  T3 - T1   0.04862 0.302 24   0.161  0.9985
##  T4 - T2   0.32214 0.302 24   1.066  0.7130
##  T4 - T1   0.26695 0.302 24   0.883  0.8134
##  T4 - T3   0.21833 0.302 24   0.722  0.8872
## 
## Collection = C2:
##  contrast estimate    SE df t.ratio p.value
##  T4 - T1   0.11588 0.302 24   0.383  0.9804
##  T2 - T1   0.14889 0.302 24   0.493  0.9600
##  T2 - T4   0.03301 0.302 24   0.109  0.9995
##  T3 - T1   0.15471 0.302 24   0.512  0.9554
##  T3 - T4   0.03883 0.302 24   0.128  0.9992
##  T3 - T2   0.00582 0.302 24   0.019  1.0000
## 
## Collection = C3:
##  contrast estimate    SE df t.ratio p.value
##  T1 - T4   0.25054 0.302 24   0.829  0.8402
##  T3 - T4   0.60346 0.302 24   1.997  0.2170
##  T3 - T1   0.35292 0.302 24   1.168  0.6523
##  T2 - T4   0.93851 0.302 24   3.106  0.0233
##  T2 - T1   0.68797 0.302 24   2.277  0.1319
##  T2 - T3   0.33504 0.302 24   1.109  0.6878
## 
## P value adjustment: tukey method for comparing a family of 4 estimates
```

```r
#filter significant comparisons
lsd_df <- as.data.frame(Results_lsmeans$comparisons) # Table form  
signif_lsmean <- lsd_df[lsd_df$p.value< 0.05,] # only  Collection 3 T2 - T4   p = 0.0233 is Significant
signif_lsmean
```

```
##  contrast Collection  estimate        SE df t.ratio p.value
##  T2 - T4  C3         0.9385082 0.3021985 24   3.106  0.0233
## 
## P value adjustment: tukey method for comparing a family of 4 estimates
```
### Boxplots to summarize Alpha diversity (mean _shannon) based on global fungal community


```r
# Create the box plot Treatments
ggplot(alpha_div_summary, aes(x=Treatment, y=mean_shannon, fill=Collection)) + 
  geom_boxplot() +
  scale_fill_manual(values=c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                    labels=c("Vegetative", "Flowering", "Seed development")) +
  scale_x_discrete(labels=c("Conventional", "No till", "Low Input", "Organic")) +
  xlab("Management") + ylab("Shannon diversity index") +
  ggtitle("Alpha Diversity Boxplot")+
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) 
```

![](Wheat_microbiomeReproScripts_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
# Potential code for adding letters to plot geom_text(data=as.data.frame(signif_lsmean), aes(x = Contrast, y = emmean, label=letters))
```
Alpha diversity conclusion
The Boxplot shows fungal community diversity peaks during flowering across all the treatments. It is lowest during vegetative period except for the Organic management where the lowest diversity was recorded during Seed development. 
The analysis of variance (ANOVA) shows that the variation in diversity is explained by the time of collection since it has a significant p-value (0.002699). The Tukey test is then used to show that during Seed development the means of No till is significantly different from the Organic management.

This analysis therefore supports our prediction that No till management supports more fungal diversity that other management strategies. Specifically during seed development, it supports more diversity than organic management


## 3.2 Beta diversity for global community to test if microbial communities were dissimilar within management types across time

```r
# Beta diversity 
fungi.dist.bray = phyloseq::distance(fungi.edited, "bray") # create bray-curtis distance matrix
fungi.ord <- ordinate(fungi.edited, "PCoA", "bray")
global.nmds <- plot_ordination(fungi.edited, ordination = fungi.ord, type = "samples") 
global.nmds.data <- global.nmds$data

# Permanova of communities considering all factors and their interactions
permanova <- adonis2(fungi.dist.bray~Collection*as.factor(Tissue)*as.factor(Treatment), as(sample_data(fungi.edited), "data.frame"), permutations = 9999) 
print(permanova) #expect PERMANOVA to show that within each management strategy, fungal communities from different growth stages (P < 0.001) had significantly different centroids.
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 9999
## 
## adonis2(formula = fungi.dist.bray ~ Collection * as.factor(Tissue) * as.factor(Treatment), data = as(sample_data(fungi.edited), "data.frame"), permutations = 9999)
##                                                    Df SumOfSqs      R2       F
## Collection                                          2    7.249 0.11743 15.7995
## as.factor(Tissue)                                   2    3.440 0.05573  7.4984
## as.factor(Treatment)                                3    2.126 0.03444  3.0891
## Collection:as.factor(Tissue)                        4    2.661 0.04310  2.8993
## Collection:as.factor(Treatment)                     6    2.585 0.04187  1.8778
## as.factor(Tissue):as.factor(Treatment)              6    1.466 0.02375  1.0651
## Collection:as.factor(Tissue):as.factor(Treatment)  12    2.514 0.04073  0.9132
## Residual                                          173   39.689 0.64294        
## Total                                             208   61.730 1.00000        
##                                                   Pr(>F)    
## Collection                                        0.0001 ***
## as.factor(Tissue)                                 0.0001 ***
## as.factor(Treatment)                              0.0001 ***
## Collection:as.factor(Tissue)                      0.0001 ***
## Collection:as.factor(Treatment)                   0.0002 ***
## as.factor(Tissue):as.factor(Treatment)            0.2993    
## Collection:as.factor(Tissue):as.factor(Treatment) 0.7895    
## Residual                                                    
## Total                                                       
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#beta diversity plot: Principal Component chart
ggplot() + 
  geom_point(data = global.nmds.data, aes(x = Axis.1, y = Axis.2, shape = as.factor(Collection), fill = as.factor(Tissue)), alpha = 0.8, size = 2) +
  theme_bw() +
  ylab("PCoA2") + 
  xlab("PCoA1") +
  scale_fill_manual(values=cbbPalette) +
  stat_ellipse(data = global.nmds.data, aes(x = Axis.1, y = Axis.2, group = Treatment), type = "norm", linetype = 2) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  facet_wrap(~ Treatment, nrow = 1, labeller = labeller(Treatment = c(T1 = "Conventional", `T2` = "No till", `T3` = "Low Input", T4 = "Organic"))) +
  ggtitle("Effect of management strategies on beta diversity")
```

![](Wheat_microbiomeReproScripts_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Beta diversity conclusion
PERMANOVA tests showed that within each management strategy, fungal communities were dissimilar across different growth stages (P < 0.001).
The findings support our prediction that there would be community shifts as the season progressed.





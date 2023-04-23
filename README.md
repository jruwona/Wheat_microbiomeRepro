# Wheat_microbiomeRepro

The microbiome dataset utilized in this project is used with permission from [Gdanetz, K., & Trail, F. 2017](https://apsjournals.apsnet.org/doi/abs/10.1094/PBIOMES-05-17-0023-R). Only the fungal ITS amplicon sequencing dataset is used for this class project.

Wheat plants were sampled from a wheat/maize/soybean crop rotation site that implements four different crop management strategies (conventional/no till/low input/organic). The plants were sampled (leaves, stems and roots) three times during the growing season and the fungal communities are analyzed. The experimental design is randomized treatments with 4 treatments (crop management) and 6 replications. 

This work aims to:
1. Determine if the global wheat fungal endophytic communities differ by management strategies across the growing period.
2. Examine the effect of management strategies on beta diversity fungal communities originating from plant organs.

Hypotheses:
1. We hypothesized that no till management supported more diverse fungal communities compared to other management options.
2. We also expected that fungal communities differ across the growth stages.


Within-sample diversity, alpha diversity statistics were calculated with the ‘phyloseq’ package, significance was tested by analysis of variance (ANOVA) and Tukey’s honest significant difference. Sample dissimilarity, beta diversity were calculated using a Bray-curtis distance matrix (Principal Coordinate Analysis). All graphs were generated with ‘ggplot2’. 

This workflow is fully reproducible by running the following scripts provided all nessesary software is installed: 
1. [Datasets](https://github.com/jruwona/Wheat_microbiomeRepro/tree/main/Data_files) 
2. [Analysis](https://github.com/jruwona/Wheat_microbiomeRepro/blob/main/Wheat_microbiomeReproScripts.Rmd)
3. [Workflow Illustration](https://github.com/jruwona/Wheat_microbiomeRepro/blob/main/Wheat_microbiomeReproScripts.md)


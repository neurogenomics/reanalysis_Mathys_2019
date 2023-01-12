# reanalysis_Mathys_2019
A re-analysis of the [Single-cell transcriptomic analysis of Alzheimer’s disease](https://www.nature.com/articles/s41586-019-1195-2) using a standardised
data processing and pseudobulk differential expression approach.

## Standardised processing protocol
We used [scFlow](https://www.biorxiv.org/content/10.1101/2021.08.16.456499v2) 
(v0.6.1) for all of the 24 patients with Alzheimer's disease pathology and the 
24 control patients. The config file with all the parameters is available in
`./scFlow_files`. This approach resulted in more stringent quality control, 
leading to the exclusion of more, problematic samples.

## Download processed data
The single-cell object (SCE) of the 
[Mathys et al.](https://doi.org/10.1038/s41586-019-1195-2) study into 
Single-nucleus transcriptomic analysis and differential expression (DE) of 
Alzheimer’s disease data after processing with scFlow is available for 
download:

```
#TODO - add wget when live
```

This will be needed to run the analysis below.

## Run analysis
The `run_reanalysis_Mathys_19.R` in the `R` folder can be used to derive the 
true DEGs from the reprocessed Mathys et al., 2019 Alzheimer's disease patient 
snRNA-Seq data. This script uses a custom written function to apply pseudobulk
differential analysis to any single-cell dataset (see `sc_cell_type_de.R`).

## Comparison against original findings
See `Mathys_results_comparison_pb.html` for a comparison between the results 
reported in the original publiciation and our findings using pseudobulk 
differential expression and a standardised data processing approach. Note the 
authors took cells as independent replicates, a cell-level analysis, in their 
work and compared the consistency in directionality and rank of their 
differentially expressed genes (DEGs) against a poisson mixed model. We show 
that these DEGs are just an artefact of taking cells as independent replicates 
by plotting the number of DEGs found against the cell counts. There is a near 
perfect correlation for them but not for us when we use pseudobulk.
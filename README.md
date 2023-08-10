# Re-analysis Mathys et al., 2019
A re-analysis of the [Single-cell transcriptomic analysis of Alzheimer’s disease](https://www.nature.com/articles/s41586-019-1195-2) using a standardised
data processing and pseudobulk differential expression approach.

## Standardised processing protocol
We used [scFlow](https://www.biorxiv.org/content/10.1101/2021.08.16.456499v2) 
(v0.6.1) for all of the 24 patients with Alzheimer's disease pathology and the 
24 control patients. The config file with all the parameters is available in
`./scFlow_files`. This approach resulted in more stringent quality control, 
leading to the exclusion of more, low quality cells.

## Download processed data
An single-cell object (SCE) of the 
[Mathys et al.](https://doi.org/10.1038/s41586-019-1195-2) study into 
Single-nucleus transcriptomic analysis and differential expression (DE) of 
Alzheimer’s disease data after processing with scFlow is available for 
download from Figshare:

```
wget https://figshare.com/ndownloader/files/38819949 -O ./data/sce.qs
```

Note this includes the processed count matrix and associated metadata. This 
will be needed to run the analysis below. Also note that here we match the cell
types to the original paper but we have also made available a more granular cell
type groupings which you can access in the metadata - the `cluster_celltype` 
column rather than the `allan_celltype` column.

## Run analysis
The `run_reanalysis_Mathys_19.R` in the `R` folder can be used to derive the 
EGs from the reprocessed Mathys et al., 2019 Alzheimer's disease patient 
snRNA-Seq data. This script uses a custom written function to apply pseudobulk
differential analysis to any single-cell dataset (see `sc_cell_type_de.R`).

## Docker file
We also provide a docker file to create an image to rerun the analysis - 
removing the need
of the user to install all the dependencies themselves. Simply run the 
following with docker installed:

```
docker pull neurogenomicslab/reanalysis_mathys_2019
```

or recreate the docker image with:

```
docker build -t reanalysis_mathys_2019 .
```

Whether you pull or recreate the image, next run it:

```
docker run -e PASSWORD=reanalysis --rm -p 8787:8787 reanalysis_mathys_2019
```

Now navigate to `localhost:8787` in a web browser and log in with:
username: rstudio
password: reanalysis
to access the docker image. Then clone the repo in the terminal of Rstudio
with:

```
git clone https://github.com/neurogenomics/reanalysis_Mathys_2019
#install data
cd reanalysis_Mathys_2019/
wget https://figshare.com/ndownloader/files/38819949 -O ./data/sce.qs
#rerun analysis
cd R/
Rscript run_reanalysis_Mathys_19.R
```

Once this runs you can look in the results folder or also knit the Rmd 
file (Mathys_results_comparison_pb.html) to view the results.

## Comparison against original findings
[See results](https://neurogenomics.github.io/reanalysis_Mathys_2019/Mathys_results_comparison_pb.html) 
### details
`Mathys_results_comparison_pb.html` gives a comparison between the results 
reported in the original publication and our findings using pseudobulk 
differential expression and a standardised data processing approach. Note the 
authors took cells as independent replicates, a cell-level analysis, in their 
work and compared the consistency in directionality and rank of their 
differentially expressed genes (DEGs) against a poisson mixed model. We show 
that these DEGs are just an artefact of taking cells as independent replicates 
by plotting the number of DEGs found against the cell counts. There is a strong 
correlation for their results but not for the pseudobulk DEGs.

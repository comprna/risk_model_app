# B-ALL risk model app

Available from [this link](https://aclosa.shinyapps.io/bALL_predictor_model/)

This interactive web application allows to calculate risk scores for patients according to the model from our [pre-print](https://www.biorxiv.org/content/10.1101/2021.12.13.472370v3).

## Input options

### Data pre-processing

Expected data is a file of normalised gene logCPMs with the gene names as rownames or as the first column of the dataframe.

Gene names and ENSEMBL gene ids are supported.

The model supports GRCh38 annotation only.

### Data format

Please indicate the following in the upload panel:
* Header: wether the data contains a header. Headers with the sample identifiers are recommended.
* Genes as rownames: if genes are the rownames of the provided dataframe (gene column does not have a column name) please tick this box. Otherwise, if the column containing gene indentifiers has a column name (i.e. Genes, GeneID, etc.) please untick this box.
* Separator: indicate the separator in your input file

Sample files for the supported input formats can be found in the sample_input_files folder. 
                    
### Missing genes

Genes that cannot be mapped will be assigned a missing value. This can affect the performance of the model. We recomend checking gene identifiers

### Score threshold
Threshold to be used to classify between high and low risk. By default 0.7 is selected (for further details see publication). Patients above the threshold should be classified as high-risk and patients below the threshold as low-risk.

## Outputs

### Gene expression

Gene expression tab contains logCPM gene expression in the provided data.

In the visualization a dashed line is added for the mean logCPM expression in the TARGET cohort, used to train the model.

### Model scores

Scores according to our risk model are displayed for the provided data. Below the Kaplan-Meyer and score distribution for the TARGET cohort are displayed according to the selected threshold.

The model is described in the publication [REF coming soon]

### Explore TARGET

Visualization of score values and expression of genes in the model according to available clinical variables for the TARGET cohort.

TARGET samples are available from  the [TARGET data portal at National Cancer institute (NIH)](https://ocg.cancer.gov/programs/target/data-matrix)
                    

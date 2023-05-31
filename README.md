# TFG: map of the repository
This is a repository of the code I wrote and the data I obtained during my Bachelor's Degree Final Project. The following file contains a map of the repository and its contents. 
##### Disclaimer: all scripts are commented for proper comprehension. However, they might contain paths relative to my local folder organization, so they may need to be adapted for future uses. 
## Data
This folder contains all the data used and generated during the project. 
### InputData
This folder contains the original files from which the candidate regions where extracted, to showcase the heterogeneity of formats they were published in and ensure accuracy. Sources are numbered from 1-12 in no particular order, but this notation is kept throughout the analysis.
### StandarizedData
This folder contains an individual CSV file for each source, corresponding to its candidate regions after standarization.
### Results
This folder contains 2 files: 
#### - AllData.csv
This is merely a merged version of all the standarized data from all sources, with the addition of gene annotation. 
#### - GenesData.csv
This file consists of a table recording each individual gene's appearence in candidate regions, for the gene-based analysis performed.
### HumanGenome.txt
A file containing a chromosome-based description of the human genome (hg38), extracted from https://www.ncbi.nlm.nih.gov/genome/?term=txid9606%5borgn%5d.

# READ ME

The order of the codes is as follows:

## 1. MIC extraction and FC construction

To extract MIC and construct functional connecome fingerprinting from a given dataset.

This step contains two R codes. *MIC_extraction.R* and *FC_construction.R* 

HCP-YA minimal pre-processed dataset is used here as an example. To analyze the other dataset, please change the file path approporiately.

## 2. Identification matrix

To construct the identification matrix by using the results from step 1

This step contains a R code, *identification_matrix.R*

MIC from HCP-YA minimal pre-processed data is used here as an example. To analyze other datasets, the tasks and file path should be replaced appropriately

## 3. Calculating the identification rate 

To calculate the identification rate from a given identification matrix (step 2)

This step contains a R code, *IDrate.R*

MIC identification matrix from HCP-YA minimal pre-processed data is used here as an example; To analyze other datasets, the tasks and file path should be replaced appropriately

## 4. Plot the brain map of MIC

This step contains a R code, *brain.R*

## 5. Sliding threshold analysis for HCP-YA ICA-FIX resting data

This step contains a R code, *threshold_analysis_HCPYA_ICAFIX.R*

## 6. histogram of explained variance for ICA-FIX data

This step contains a R code, *Varhist.R*

## 7. thresholding for HCP-YA minimal pre-processed data

This step contains a R code, *threshold_analysis_HCPYAminimal.R*

threshold of 0.65 is used here as an example

## 8. denoising pipeline

this step contains a python code (UNIX virtual system), *denoise_pipeline.ipynb*
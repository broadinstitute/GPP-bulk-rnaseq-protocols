# RNAseq Alignment and Quantification

See [`Method`](#Method) section for a detailed explanation of how alignment and quantification is carried out.    
See [`Protocol`](#Protocol) section for instructions.    


# Method

# Protocol
## Requirements
* Access to [Terra](https://terra.bio/)
* Google Cloud CLI

### Install Google Cloud CLI with Homebrew on MacOS
* Install x-code if you don't already have it 
	* `xcode-select --install`
* Install [Homebrew] (https://brew.sh/) as the MacOS package manner if you don't already have it
* Install the CLI
	* `brew install --cask google-cloud-sdk`
* Log into your Google Cloud account
	* `gcloud auth login`

### Create a Terra sample file
1. Make a *tab delimited* sample file with the headers : `entity:sample_id`, `fastq1`, `fastq2`.   
2. Fill the `fastq1` column with the R1 files and `fastq2` with R2 (While you can do this manually by clicking through every file, we recommend using the Google Cloud CLI to list the files). 
	* Find the Google bucket path to where your data files are located (e.g. `gs://gpp_rnd_rnaseq/TP53`).   
	* Get R1 reads `gcloud storage ls gs://path/to/folder | egrep '\.gz$' | grep R1`
	* Get R2 reads `gcloud storage ls gs://path/to/folder | egrep '\.gz$' | grep R2`
3. Determine a unique and readable sample name for the sample ID column.   
	* Helpful to prefix all your samples with a unique identifier so you can search for them easily.    
4. Save the file with`.tsv `or a `.txt` extension. Terra will accept either one.


### Add samples to Terra
1. Click on the **Data** tab of the broad-gpp/RNASeq_analysis workspace.      
2. Click on **Import Data** -> Upload TSV.  
![Import Data Screenshot](https://user-images.githubusercontent.com/7750862/217667602-da01c04a-9d16-42a8-adbe-a32746b56e39.png) 
3. Upload the TSV/TXT file you created in the last step and click **"Start Import Job"**.  
4. Make sure the number of samples added matches the number in your file.   

### Start a new job
1. Click on the **Workflows** tab of the broad-gpp/RNASeq_analysis workspace.   
2. Click on **RNA_pipeline**.   
3. Click on **"Select Data"** in Step 2 and select your samples. *Leave all the default settings as-is*.
	![Screenshot 2023-02-08 at 5 36 19 PM](https://user-images.githubusercontent.com/7750862/217667886-3e382e85-f4ab-493f-a5ee-94e457cff8d6.png)
4. Click **"Run Analysis"**
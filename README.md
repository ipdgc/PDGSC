# PDGSC Analysis Pipeline

## Table of contents
- [Introduction](#introduction)
- [Summary of the PDGSC pipeline](#pipeline)
- [Useful scripts for exome analyses](#scripts)
- [Some general notes for beginners running analyses on the cloud](#cloud)


## Introduction <a name="introduction"></a>

Collection of scripts used in the analysis and QC of the PDGSC data, as well as some useful scripts for other analyses.

The PDGSC is funded by MJFF => https://www.michaeljfox.org/grant/parkinsons-disease-genetics-sequencing-consortium


## Summary of the PDGSC pipeline <a name="pipeline"></a>

To reproduce the PDGSC analysis pipeline, do the following:

1. Clone this github repository

2. Install required/useful tools
````bash
bash ./pdgsc/scripts/shell/install_tools.sh
````

3. Download analysis files (VCF etc)
````bash
bash ./pdgsc/scripts/shell/download_analysis_files.sh
````

4. Generate depth metrics and calculate individual-level depth statistics
````bash
bash ./pdgsc/scripts/shell/get_ad_info.sh
````

5. Filter out poor quality samples and perform GWAS-style individual-level QC
````bash
bash ./pdgsc/scripts/shell/individual_qc.sh
````

6. Match ADSP controls to samples from cohorts with no controls
````bash
bash ./pdgsc/scripts/shell/casecontrol_matching.sh
````

7. Perform variant QC and generate covariates for analysis
````bash
bash ./pdgsc/scripts/shell/variant_qc_generate_covariates.sh
````

8. Run Rvtests to generate covariance files
````bash
bash ./pdgsc/scripts/shell/run_rvtests.sh
````

9. Run Raremetals to meta-analyse
````bash
bash ./pdgsc/scripts/shell/run_raremetals.sh
````

10. Do some post-processing
````bash
bash ./pdgsc/scripts/shell/subset_results_postqc.sh
````

That's it!

Flowchart etc here

## Useful scripts for exome analyses <a name="scripts"></a>

If you're looking to perform an analysis on a candidate gene(s) of interest, then the following script might be useful. You can use it to extract the gene of interest, keep only samples and variants that are covered well in the gene, generate an annotated table of all the variants found in the gene, and perform burden tests on the gene.

````bash
bash ./pdgsc/scripts/shell/extract_annotated_variants_for_gene_of_interest.sh
````

## Some general notes for beginners running analyses on the cloud <a name="cloud"></a>

To connect to the Google Cloud using your computer, first install the Google Cloud SDK on your system by following the instructions for your operating system here:
https://cloud.google.com/sdk/docs/downloads-interactive

Then, to log in and initialise your account, run:

````bash
gcloud init
````

You should receive a prompt to log in to your Google account, and to connect your account to the PDGSC project. If this doesn't happen, you might need to ask to be added to the project.

Once you are all set up, you can start a virtual machine with the following command (just pick a name for your machine, how much hard disk size you think you'll need, and a machine type):

A list of machine types and their specs can be found here:
https://cloud.google.com/compute/docs/machine-types

````bash
VM_NAME="type your machine name here"
DISK_SIZE="type your disk size here"
MACHINE_TYPE="type your machine type here"

gcloud compute instances create ${VM_NAME} --zone us-central1-f --image-family ubuntu-1804-lts --image-project ubuntu-os-cloud  --machine-type ${MACHINE_TYPE} --maintenance-policy MIGRATE --boot-disk-size ${DISK_SIZE} --boot-disk-type pd-standard --boot-disk-device-name ${VM_NAME}
````
Or if you are using the Google SDK on Windows, then:
````dos
SET VM_NAME="type your machine name here"
SET DISK_SIZE="type your disk size here"
SET MACHINE_TYPE="type your machine type here"

gcloud compute instances create %VM_NAME% --zone us-central1-f --image-family ubuntu-1804-lts --image-project ubuntu-os-cloud  --machine-type %MACHINE_TYPE% --maintenance-policy MIGRATE --boot-disk-size %DISK_SIZE% --boot-disk-type pd-standard --boot-disk-device-name %VM_NAME%
````

Once you've successfully set up a virtual machine, you can connect to it with by running:

````bash
gcloud compute ssh ${VM_NAME}
````

Or on windows:

````dos
gcloud compute ssh %VM_NAME%
````

And you're good to go! Once you are done with your analysis, remember to delete your virtual machine either by navigating to the virtual machines section on your Google Cloud Dashboard, or by running the following inside the virtual machine:

````bash
gcloud compute instances delete ${VM_NAME}
````

It might be a good idea to add this to the end of your analysis script, so that the virtual machine is deleted once your analysis is finished and is not left running for no reason.
Remember to save the output/results of your analysis before you do this though! You can do this by uploading the files into the bucket by running:

````bash
OUTPUT_TO_BE_SAVED="name of the file you want to save in the bucket"
BUCKET_ADDRESS="the location in the bucket where you want to save your file"

gsutil -m cp ${OUTPUT_TO_BE_SAVED} ${BUCKET_ADDRESS}
````

If you're planning to run scripts from this repository, it might be easier to download this repository to the VM, by running the following. It might also be useful to install some basic useful tools that are often needed for analysis after setting up the VM (you can edit the install_tools.sh script to include/remove tools according to your needs):

````bash
git clone https://github.com/ipdgc/pdgsc.git
bash ./pdgsc/scripts/shell/install_tools.sh
````

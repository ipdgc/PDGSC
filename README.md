# PDGSC Analysis Pipeline

## Introduction

Collection of scripts used in the analysis and QC of the PDGSC data, as well as some useful scripts for other analyses.

## Summary of the PDGSC pipeline

Flowchart etc her

## Useful scripts for exome analyses

Some useful scripts for secondary analyses

## Some general notes for beginners running analyses on the cloud:

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

Once you've successfully set up a virtual machine, you can connect to it with by running:

````bash
gcloud compute ssh ${VM_NAME}
````

And you're good to go! Once you are done with your analysis, remember to delete your virtual machine either by navigating to the virtual machines section on your Google Cloud Dashboard, or by running:

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

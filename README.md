# PDGSC Analysis Pipeline

## Table of contents
- [Introduction](#introduction)
- [Summary of the PDGSC pipeline](#pipeline)
- [Useful scripts for exome analyses](#scripts)
- [Some general notes for beginners running analyses on the cloud](#cloud)
- [PDGSC Acknowledgement](#PDGSC)


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



## PDGSC Acknowledgement <a name="PDGSC"></a>



PDGSC (Parkinson Disease Genetic Sequencing Consortium) is a collaborative group of investigators working in the area of PD genetics through the analysis of high content sequencing. The PDGSC has been supported by The Michael J. Fox Foundation for Parkinson’s Research, the National Institute on Aging Intramural Research Program, and the National Institute of Neurological Disorders and Stroke. A full list of the participants support is provided below.

The members of the consortium are: Marco Abreu (Indiana University School of Medicine, USA), Gary W. Beecham (John P. Hussman Institute for Human Genomics, University of Miami, USA), Sara Bandres-Ciga (National Institute on Aging, USA), Cornelis Blauwendraat (National Institute on Aging and National Institute of Neurological Disorders and Stroke, USA), Jose Bras (University College London, UK), Alexis Brice (Brain and Spine Institute (ICM), France), Kathrin Brockmann (University of Tübingen/DZNE, Germany), Zeynep Akdemir (Baylor College of Medicine, USA), Patrick F Chinnery (University of Cambridge, UK), Jean-Christophe Corvol (Brain and Spine Institute (ICM), France), Fabrice Danjou (Brain and Spine Institute (ICM), France), Aaron Day-Williams (MRL, Merck & Co., Inc., Boston, MA, USA), John D Eicher (MRL, Merck & Co., Inc., Boston, MA, USA), Karol Estrada (Biogen, USA), Daniel M Evans, Faraz Faghri (National Institute on Aging, USA, University of Illinois at Urbana-Champaign, USA), Samuel Evetts (University of Oxford, UK), Ilaria Guella and Matthew J Farrer (Centre for Applied Neurogenetics, University of British Columbia, Canada), Tatiana Foroud (Indiana University School of Medicine, USA), Steve Finkbeiner (Gladstone Institutes/UCSF, USA), Thomas Gasser (University of Tübingen/DZNE, Germany), J Raphael Gibbs (National Institute on Aging, USA), John Hardy (University College London, UK), MTM Hu (University of Oxford, UK), Joseph Jankovic (Baylor College of Medicine, USA), Hallgeir Jonvik (University College London, UK), Demis A Kia (University College London, UK), Christine Klein (Institute of Neurogenetics, University of Luebeck, Germany), Rejko Krüger (Luxembourg Centre for Systems Biomedicine, Luxembourg), Dongbing Lai (Indiana University School of Medicine, USA), Suzanne Lesage (Brain and Spine Institute (ICM), France), Christina M. Lill (Institute of Neurogenetics, University of Luebeck, Germany), Steven J. Lubbe (Ken and Ruth Davee Department of Neurology, Northwestern University, Chicago, IL, USA), Timothy Lynch (Dublin Neurological Institute, Mater Misericordiae University Hospital, Ireland), Kari Majamaa (University of Oulu, Finland), Eden R. Martin (John P. Hussman Institute for Human Genomics, University of Miami, USA), Patrick May (Luxembourg Centre for Systems Biomedicine, Luxembourg), Brit Mollenhauer (University Medical Center Goettingen, Department of Neurology, Germany), David Murphy (University College London, UK), Huw R Morris (University College London, UK), Mike A Nalls (National Institute on Aging, USA, Data Tecnica International, USA), Khanh-Dung Nguyen (Biogen, USA), Karen Nuytemans (John P. Hussman Institute for Human Genomics, University of Miami, USA), Lasse Pihlstrom (Oslo University Hospital, Norway), Alan Pittman (University College London, UK), Lea R'Bibo (University College London, UK), Laurie Robak (Baylor College of Medicine, USA), Owen A. Ross (Mayo Clinic Jacksonville, USA), Cynthia Sandor (University of Oxford, UK), Barbara Schormair (Institute of Neurogenomics, Helmholtz Zentrum München - Deutsches Forschungszentrum für Gesundheit und Umwelt (GmbH), Munich, Germany), William K. Scott (John P. Hussman Institute for Human Genomics, University of Miami, USA), Manu Sharma (Centre for Genetic Epidemiology, Institute for Clinical Epidemiology,University of Tubingen, Germany), Joshua M Shulman (Baylor College of Medicine, USA), Ari Siitonen (University of Oulu, Finland), Javier Simón-Sánchez (University of Tübingen/DZNE, Germany), Andrew B Singleton (National Institute on Aging, USA), David J Stone (MRL, Merck & Co., Inc., West Point, PA, USA), Konrad Szewczyk-Krolikowski (University of Oxford, UK), Manuela MX Tan (University College London, UK), Paul Tomlinson  (University of Oxford, UK), Mathias Toft (University of Oslo, Norway), Richard Wade-Martins (University of Oxford, UK), Claudia Trenkwalder (Dept. Neurosurgery, University Medical Center, Goettingen, Germany), Caleb Webber (University of Oxford, UK), Wei Wei (University of Cambridge, UK), Jeffery M. Vance (John P. Hussman Institute for Human Genomics, University of Miami, USA), Nigel M Williams (Cardiff University, UK), Juliane Winkelmann (Institute of Neurogenomics, Helmholtz Zentrum München - Deutsches Forschungszentrum für Gesundheit und Umwelt (GmbH), Munich, Germany), Zbigniew K. Wszolek (Mayo Clinic Jacksonville, USA), Pauli Ylikotila (University of Turku, Finland), Alexander Zimprich (Department of Neurology, Austria).

For correspondence regarding PDGSC please contact: Andrew B Singleton, Laboratory of Neurogenetics, National Institute on Aging, National Institutes of Health, Bethesda, USA; singleta@mail.nih.gov

The work of the PDGSC was supported by: The Intramural Research Program of the National Institute on Aging, National Institutes of Health, part of the Department of Health and Human Services (ZO1 AG000949), The Extramural Research Program of the National Institutes of Health (R01NS096740, R01NS037167, R01NS078086, P50NS071674, P50NS072187), The Michael J. Fox Foundation for Parkinson’s Research, The Department of Defense (USAMRAA, the Mitochondrial Disease (W81XWH-17-1-0249) and the Parkinson's Research Program managed through CDMRP), the National Institute of Neurological Disorders and Stroke, the Canadian Consortium for Neurodegeneration in Aging, the Canada Excellence Research Chairs program, the Canadian Institutes of Health Research,  the Hermann and Lilly Schilling Foundation, the German Research Foundation (FOR2488), the German Federal Ministry of Education and Research (BMBF) under the funding code 031A430A and the EU Joint Programme -Neurodegenerative Diseases Research (JPND) project under the aegis of JPND (www.jpnd.eu) through Germany, BMBF, funding code 01ED1406, Sigrid Juselius Foundation, the Medical Research Council UK (MC_UP_1501/2 & 13044), Parkinson’s UK (grants 8047, J-0804, G-1502), the Medical Research Council UK (G0700943, G1100643) the Wellcome Trust (101876/Z/13/Z), the Fonds National de Recherche (FNR; NCER-PD) Luxembourg.



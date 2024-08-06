# Base Editing Sensor Analysis Pipeline
Author: Sam Gould | Last Updated: 8/6/24

This repository provides a step-by-step breakdown of going from receiving your sequencing data to getting a processed sensor dataset that includes (1) sequencing quality breakdown, (2) guide counts, (3) sensor editing analysis.

This pipeline is based on the following sequencing strategy:
![Sensor](images/seq_strategy.png)

**This pipeline would need to modified to work with alternative sequencing strategies.**

## Downloading your sequencing data
- Log into your Luria cluster account by going into the Terminal and entering your login credentials. This follows the format:

```  
Ssh username@luria.mit.edu
* click enter *
Password
* click enter *
```  
- To download the data to *your server folder*, use the instructions provided in the email from the sequencing core. For example: 


"1. For users with BioMicro Linux cluster (luria.mit.edu) accounts only:
Run rsync in srun or a sbatch script to sync the data to your own data folder (e.g. ~/data/my_target_dir)"
```
srun rsync -av /net/bmc-pub17/data/bmc/public/Sanchez/project_name ~/data/my_target_dir”
```
- You need to replace the folder where the data will be stored, e.g.: 
```
srun rsync -av /net/bmc-pub17/data/bmc/public/Sanchez/project_id /net/bmc-lab2/data/lab/sanchezrivera/your_folder_name
```
*your_folder_name* = the name of your existing folder on the cluster; this step won't work otherwise. For me this is "samgould".


## Step 0: Creating conda environments

In order to run these scripts, you need to create **2** conda environments. Conda environments are essentially sandboxes that allow python scripts to run while referencing all of their required packages/package versions. These environments allow the scripts to run on the cluster, otherwise you would get errors of "package not installed" when trying to import the packages.

The .yml files for creating these conda environments are available in the **conda_environments** folder of this repository. 

To create these conda environments, **copy these .yml files to your own folder on the server**, and then run the following commands after logging into the cluster:

```   
#first start an interactive session outside the head node
srun --pty bash

#cd into your server folder where .yml files are (replace my name)
cd /net/bmc-lab2/data/lab/sanchezrivera/samgould

#load conda
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit

#and then create the two enviornments
#each of these steps may take a while

conda env create -n sensor_lib_sg -f ./sensor_lib_sg.yml
conda env create -n crispresso_env -f ./crispresso_environment.yml


#you can validate these are installed with this command:
conda env list
```  

More info is available about troubleshooting conda environments at the [Luria website](https://igb.mit.edu/mini-courses/advanced-utilization-of-igb-computational-resources/package-management/conda-environments#sharing-conda-environments).


## Step 0.1: 



| Syntax | Description |
| --- | ----------- |

- “gRNA_id”: a unique identifier for each guide RNA
- “Protospacer”: G+19 protospacer, this includes the “G” start.
- **Hamming_BC**: the 15 nt unique sensor barcode
- **sensor_wt**: wildtype 42 nt sensor sequence
- **sensor_alt**: correctly edited 42 nt sensor 


![Step 1](images/1.png)
![Step 2](images/2.png)
![Step 3](images/3.png)

## Post-processing

## Addendum: MAGeCK

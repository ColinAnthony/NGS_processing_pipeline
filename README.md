# NGS_processing_pipeline

## Overview:
### This pipeline
1) Creates the required file structure
2) Demultiplexes sequences based on the presence of a primer
3) Does other things

### Setup process:

 1) Clone the git repository
 2) etc
 

### Dependecies:
Data from PrimerID workflow

    1) Motifbinner
    2) contam repo
    3) pyhton3
    4) scipy?
    5) pandas
    6) numpy
    7) Biopython?  
  

### To run:

```python3 demultiplex.py --config_file <config.json> --output_dir <output directory>```



### Using the config file:

The configuration for the run is stored in a JSON file for reproducibility. An example file is provided 
which can be edited to suit your run. 

**The JSON file contains:** 

  **input_data**

* fwd_fastq_file - This is the fastq file containing the forward (R1) sequences.
* rev_fastq_file - This is the fastq file containing the reverse (R2) sequences.
* primer_csv - This file contains the primers and their associated gene regions, as well as information on
whether the sequences are overlapping, and the length of sequence preceding the primer in the fastq file. 
* patient_list - This is a python list containing the patient / sample name.
* out_folder - This is the full path to the output folder.

  
  **pipelineSettings**
  
* out_prefix: The prefix name of your outfile.
* frame: The reading frame (1, 2 or 3).
* stops: Remove sequences with stop codons?
* min_read_length: The minimum read length.
* run_step: The step at which to resume the analysis if it is interrupted. 

  
  **haplotype_settings**
  
* infile: The path and name of the aligned fasta file.
* field: The field that differentiates your samples/time points (use the last field if multiple. ie: 4 for 'CAP177_2000_004wpi_V3C4_GGGACTCTAGTG_28, or 2 for SVB008_SP_GGTAGTCTAGTG_231).
* script_folder: The path to the folder containing the pipeline scripts.

 
### Output:

The output of this pipeline is X in Y format

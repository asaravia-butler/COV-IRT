## Conda Environments Containing Tools Used in the RNAseq Processing Pipelines

> **This directory holds yaml files and instructions to install conda environments containing the tools needed to process RNAseq data with the COV-IRT pipelines shown is this Repository.**

---

### Install prerequisites

  * **Anaconda**  
    To install conda environments, you'll first have to install [Anaconda](https://www.anaconda.com/). Click [here](https://docs.anaconda.com/anaconda/install/) for installation instructions.

<br>

### Install and activate COV-IRT RNAseq conda environments

  Install each conda environment by downloading the respective yaml file then runnning the following command:

  ```
  conda env create -f conda_env_name.yml
  ```

  Activate the conda environment with the following command:
   
  ```
  conda activate conda_env_name
  ``` 
  Note: Replace conda_env_name in the above commands with the name of the specific environment you wish to install. The following environments are available:
  - COVIRT_GATK
  - COVIRT_HTStream
  - COVIRT_RNAseq_R_tools
  - COVIRT_fastq_to_alignment
  - COVIRT_fgbio
  - COVIRT_gtfToGenePred
  - COVIRT_nextflow
  - COVIRT_rsem
  
  Deactivate an active conda environment with the following command:
   
  ```
  conda deactivate 
  ``` 

<br>

### Troubleshooting conda environment installation

  If you run into any issues while installing the RNAseq conda environments, update conda using the following command then try re-creating the RNAseq conda environments using the installation instructions above.
  ```
  conda update conda
  ```


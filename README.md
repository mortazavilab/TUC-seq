
## Pipeline description and installation
This repository consists of scripts that were used to process TUC-seq samples for analysis and comparison to long-TUC-seq data and for generating the graphs used in our paper.

In order to run the pipeline, you need to first place your pre-processed TUC-seq fastq files in "data/" directory and ensure that your genome reference and annotation files are deposited in the "ref/" directory.

This pipeline requires the following packages to run properly:

* [STAR](https://github.com/alexdobin/STAR)
* [eXpress](https://pachterlab.github.io/eXpress/overview.html)

Once the neccesary packages are installed and data and references are in place. First adjust the path to this repository in all the scripts and run the following scripts:

1. Run alignment.sh on the fastq file for each of the files to get the aligned files.
1. Run substitution-annotator.sh on each of the aligned sam files to categorize the reads based on the number of T to C substitutions in each read.

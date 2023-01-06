# Toxinpred2
A method for predicting toxicity of the proteins

# Introduction
ToxinPred2.0 is developed for predicting, mapping and scanning toxic peptides. More information on ToxinPred2 is available from its web server http://webs.iiitd.edu.in/raghava/toxinpred2. This page provide information about standalone version of ToxinPred2.

# Standalone

Standalone version of ToxinPred2 is written in python3 and the following libraries are necessary for a successful run:

- scikit-learn
- Pandas
- Numpy
- blastp

# Important Note

- Due to large size of the model file, we have not included it in the zipped folder or GitHub repository, thus to run standalone successfully you need to download model file and then unzip them.
- Make sure you extract the download zip files in the directory where main execution file i.e. toxinpred2.py is available.
- To download the model folder click [here].(https://webs.iiitd.edu.in/raghava/toxinpred2/stand.html)

**Minimum USAGE** 

To know about the available option for the standalone, type the following command:
```
toxinpred2.py -h
```
To run the example, type the following command:
```
toxinpred2.py -i peptide.fa

```
where peptide.fa is a input FASTA file. This will predict toxic peptides in FASTA format. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma separated variables).

**Full Usage**: 
```
Following is complete list of all options, you may get these options
usage: toxinpred2.py [-h] 
                     [-i INPUT]
                     [-o OUTPUT]
                     [-t THRESHOLD]
                     [-m {1,2}] 
                     [-d {1,2}]
```
```
Please provide following arguments

optional arguments:

  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence in FASTA format or
                        single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.6
  -m {1,2}, -- model Model
                        Model: 1: AAC based RF, 2: Hybrid, by default 1
  -d {1,2}, --display {1,2}
                        Display: 1:Toxin peptide, 2: All peptides, by
                        default 1

```

**Input File**: It allow users to provide input in two format; i) FASTA format (standard) (e.g. peptide.fa) and ii) Simple Format. In case of simple format, file should have one peptide sequence in a single line in single letter code (eg. peptide.seq). 

**Output File**: Program will save result in CSV format, in case user do not provide output file name, it will be stored in outfile.csv.

**Threshold**: User should provide threshold between 0 and 1, please note score is proportional to toxic potential of peptide.

**Models**: In this program, two models have been incorporated;  
  i) Model1 for predicting given input peptide/protein sequence as toxic and non-toxic peptide/proteins using Random Forest based on amino-acid composition of the peptide/proteins; 

  ii) Model2 for predicting given input peptide/protein sequence as toxic and non-toxic peptide/proteins using Hybrid approach, which is the ensemble of Random Forest+ BLAST+ MERCI. It combines the scores generated from machine learning (RF), MERCI, and BLAST as Hybrid Score, and the prediction is based on Hybrid Score.


ToxinPred2 Package Files
=======================
It contain following files, brief description of these files given below

INSTALLATION  	: Installation instructions

LICENSE       	: License information

envfile : This file provide the path information for BLAST and MERCI commands ,and data 
          required to run BLAST and MERCI

Database: This folder contains the blast database

progs : This folder contains the program to run MERCI

README.md     	: This file provide information about this package

toxinpred2.py 	: Main python program 

RF_model        : Model file required for running Machine-learning model

peptide.fa	: Example file contain peptide sequences in FASTA format

peptide.seq	: Example file contain peptide sequences in simple format

protein.fa	: Example file contain protein sequences in FASTA format 

# Reference
Sharma N, Naorem LD, Jain S, Raghava GPS (2022) ToxinPred2: an improved method for predicting toxicity of proteins. <a href="https://pubmed.ncbi.nlm.nih.gov/35595541/">Brief Bioinform. doi: 10.1093/bib/bbac174.</a>

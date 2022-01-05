# ToxinPred 2.0 : A webserver for designing and predicting toxic and non-toxic proteins/ peptides

# Introduction
ToxinPred 2.0 is an updated version of our server ToxinPred, developed in 2013. It is developed for predicting, mapping and scanning toxic peptides. More information on ToxinPred2.0 is available from its web server http://webs.iiitd.edu.in/raghava/toxinpred2. This page provide information about standalone version of ToxinPred2.0. Please read/cite following paper for complete information including algorithm behind ToxinPred 2.0.

**Models**: In this program, two models have been incorporated;  
i) Model1 for predicting given input peptide/protein sequence as toxic and non-toxic peptide/proteins using Random Forest based on amino-acid composition of the peptide/proteins; 

ii) Model2 for predicting given input peptide/protein sequence as toxic and non-toxic peptide/proteins using Hybrid approach, which is the ensemble of Random Forest+ BLAST+ MERCI. It combines the scores generated from machine learning (RF), MERCI, and BLAST as Hybrid Score, and the prediction is based on Hybrid Score.

**Modules/Jobs**: This program implements three modules; 
  i) Prediction: It allows users to submit multiple sequences in FASTA format, server will predict toxic and non-toxic proteins/ peptides from their primary sequence.
  ii) Motif Scan: It allows users to scan or map motifs in the query sequence using MERCI. Here, the user can search for the crucial motifs present in the query protein sequence to find out its role as toxic or non-toxic. 
  iii) BLAST Search: It allows users to search the query protein sequence against the database of known toxic and non-toxic peptides using similarity-based search method.
  
**Minimum USAGE**: Minimum usage is "toxinpred2.py -i peptide.fa" where peptide.fa is a input FASTA file. This will predict toxic peptides in FASTA format. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma separated variables).

**Full Usage**: Following is complete list of all options, you may get these options by "toxinpred2.py -h".

```

usage: toxinpred2.py [-h] -i INPUT [-o OUTPUT] [-t THRESHOLD] [-m {1,2}] [-d {1,2}]

Please provide following arguments

optional arguments:
  -h, --help            show this help message and exit
  
  -i INPUT, --input INPUT (Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code)
                        
  -o OUTPUT, --output OUTPUT (Output: File for saving results by default outfile.csv)
                        
  -t THRESHOLD, --threshold THRESHOLD (Threshold: Value between 0 to 1 by default 0.6)
                        
  -m {1,2}, -- model Model (Model: 1: AAC based RF, 2: Hybrid, by default 1)
                        
  -d {1,2}, --display {1,2} (Display: 1:Toxin peptide, 2: All peptides, by default 1)

```

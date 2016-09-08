# MultiTool - A handy set of tools for any bioinformatics-related lab!


## Table of Contents
- [Installation]()
- [Documentation]()
- [Support](#support)
 - [Help](#configuration-issueshelp)
 - [Bugs](#bugs--issues)
 - [Feature Requests](#feature-requests)
 - [Pull Requests](#pull-requests)
- [Features](#features)
- [Authors](#authors)

The project is currently setup in two main branches:
- `dev` also known as `beta` - This is where the latest features are, but you may also experience some issues with path development.
- `master` also known as `stable` - 

## Support

### Configuration issues/help
If you need any help please contact one of the authors via email.

###[Bugs / Issues]()
If you discover a bug in the script, please [create a new issue]().

###[Feature Requests]()
If you have a feature you would like added, please [create a new request]().


###[Pull Requests]()
If you'd like to make your own changes ensure your PR is made against the 'dev' branch

## Features
- [x] Gene CDS Extraction from Genbank files - Produces multifasta
- [x] Calculates gene copy numbers based on depth of coverage results from Bowtie2 alignments (Avg. Cov. & Stdev)
- [x] Find and extract sequences w/ header from a multifasta based on header patterns (e.g. Gene name)
- [x] Test reads for quality based on alignment rate, Avg. Depth of Cov & Stdev, and TS/TV. 
- [x] Calculate Proteome Size from a Proteomic or Transcriptomic Genbank file.
- [x] Count number of bases for individual fasta or multifasta.
- [x] Added a ToolBox to facilitate the lives of non-savvy terminal users (Read Documentation)
- [x] Limit the step to farm specific area for pokestops
- [x] SNP analysis based on output from Bowtie2 and Bwa-mem alignments coupled with SAMtools and GATK (w/ Piccard Tools) cross analysis.
- [ ] PSMC Population Plot Analysis


## Authors
- [Charles Sanfiorenzo Cruz]() 
- [Jenelys Ruiz Ortiz]() 


## Contributors
 
  

## Disclaimer
This program takes advantage of a lot of software tools debeloped by many academic institutions. Make sure to cite their use in accordance to their policies!

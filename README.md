# MultiTool - A handy set of tools for any bioinformatics-related lab!


## Table of Contents
- [Installation]()
- [Documentation](https://github.com/CharlesSanfiorenzo/Bioinformatics/tree/master/docs)
- [Support](#support)
 - [Help](#configuration-issueshelp)
 - [Bugs](#bugs--issues)
 - [Feature Requests](#feature-requests)
 - [Pull Requests](#pull-requests)
- [Features](#features)
- [Authors](#authors)

The project is currently setup in two main branches:
- `dev` also known as `beta` - This is where new features are tested; users may experience some issues. Dev branch repository may be found [here]()
- `master` also known as `stable`   

## Support

### Configuration issues/help
If you need any help please contact one of the [authors](#authors) via email.

###[Bugs / Issues](https://github.com/CharlesSanfiorenzo/Bioinformatics/issues)
If you discover a bug in the script, please [create a new issue](https://github.com/CharlesSanfiorenzo/Bioinformatics/issues/new).

###[Feature Requests](https://github.com/CharlesSanfiorenzo/Bioinformatics/labels/Feature%20Request)
If you have a feature you would like added, please [create a new request](https://github.com/CharlesSanfiorenzo/Bioinformatics/issues/new).


###[Pull Requests]()
If you'd like to make your own changes ensure your PR is made against the ['dev']() branch.

## Features
- [x] Gene CDS Extraction from Genbank files - Produces multifasta
- [x] Calculates gene copy numbers based on depth of coverage results from Bowtie2 alignments (Avg. Cov. & Stdev)
- [x] Find and extract sequences w/ header from a multifasta based on header patterns (e.g. Gene name)
- [x] Test reads for quality based on alignment rate, Avg. Depth of Cov & Stdev, and TS/TV. 
- [x] Calculate Proteome Size from a Proteomic or Transcriptomic Genbank file.
- [x] Count number of bases for individual fasta or multifasta.
- [x] Added a ToolBox to facilitate the lives of non-savvy terminal users (Read Documentation)
- [x] SNP analysis based on output from Bowtie2 and Bwa-mem alignments coupled with SAMtools and GATK (w/ Piccard Tools) cross analysis.
- [x] Assays Mutational Load by calculating pKa/Ks ratio using the Nei-Gojobori method.
- [ ] PSMC Population Plot Analysis
- [ ] Install batch for all required software tools.
- [ ] GUI for Non-terminal users


## Authors
- [Charles Sanfiorenzo Cruz]() - charles.sanfiorenzo@upr.edu
- [Jenelys Ruiz Ortiz]() - jenelys.ruiz@upr.edu


## Contributors
 
  

## Disclaimer
This program takes advantage of a lot of software tools developed by many academic institutions. Make sure to cite their use in accordance to their policies!

## Citation
Sanfiorenzo C. J. & Ruiz J. (2016). MultiTool - A Bioinformatic Toolset [Computer Software]. University of Puerto Rico Rio Rio Piedras. Retrieved from: https://github.com/CharlesSanfiorenzo/Bioinformatics.git

# MultiTool - A handy set of tools for any bioinformatics-related lab!


## Table of Contents
- [Installation](#installation)
- [General Usage](#general-usage)
- [Documentation](https://github.com/CharlesSanfiorenzo/Bioinformatics/MultiTool/tree/master/docs/GeneralUse.md)
- [Support](#support)
 - [Help](#configuration-issueshelp)
 - [Bugs](#bugs--issues)
 - [Feature Requests](#feature-requests)
 - [Pull Requests](#pull-requests)
- [Features](#features)
- [Requirements](#requirements)
- [Authors](#authors)

The project is currently setup in two main branches:
- `dev` also known as `beta` - This is where new features are tested; users may experience some issues. Dev branch repository may be found [here]()
- `master` also known as `stable`   

## Installation

To install, simply download all of the listed requirements and install them onto a user created directory (e.g. 'BioinformaticTools').
Then download this repository and change the path variable within 'MultiTool.sh' to the location of the directory containing the requirements. In Unix systems you could use gedit (Genome Editor) or the vi editor, but in Windows you could do this by opening 'MultiTool.sh' with the Notepad app.

## General Usage

Run MultiTool.sh and a menu will pop up:

    $ ./MultiTool.sh 
    1) Extract         4) ReadQuality    7) BaseCounter   10) ToolBox
    2) GeneCopies      5) ReadFilter     8) SNPCalling    11) Help
    3) FindGene        6) ProteomeSize   9) MutationLoad  12) Quit
    Enter your choice (use number):
    Enter your choice (use number): 8
    1) Bowtie2
    2) Bwa-mem
    Enter your choice (use number): 1
    1) SingleRead
    2) ReadPair
    Enter your choice (use number): 2
    Select your reference (.fna) and read (.fastq) files; include path.
    Note: You may choose from the following files in this directory and/or 
    subdirectories 
    #Fasta list will show here!
    Reference :
    #Fastq list will show here!
    Read_1 : 
    Read_2 : 
    And so on ... :)

For information about each individual feature, please read the documentation.

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
- [x] Filters reads using quality scores and read length 
- [x] Calculate Proteome Size from a Proteomic or Transcriptomic Genbank file.
- [x] Count number of bases for individual fasta or multifasta.
- [x] Added a ToolBox to facilitate the lives of non-savvy terminal users (Read Documentation)
- [x] SNP analysis based on output from Bowtie2 and Bwa-mem alignments coupled with SAMtools and GATK (w/ Piccard Tools) cross analysis.
- [x] Assays Mutational Load by calculating pKa/Ks ratio using the Nei-Gojobori method.
- [ ] PSMC Population Plot Analysis
- [ ] Install batch for all required software tools.
- [ ] GUI for Non-terminal users

## Requirements
* **Languages**
 * Python 2.7 or Python 3+ : https://www.python.org/
 * Perl : https://www.perl.org/get.html
 * **If on Windows**
   * Linux shell emulator : https://www.cygwin.com/
    * GNU Packages : http://gnuwin32.sourceforge.net/
* **Tools**
 * Bowtie2 : http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
 * Bwa-mem : http://bio-bwa.sourceforge.net/
 * SAMtools : http://samtools.sourceforge.net/
 * VcfTools : https://vcftools.github.io/index.html
 * Feature Extract (gb2tab) : http://www.cbs.dtu.dk/services/FeatureExtract/download.php 
 * GATK : https://software.broadinstitute.org/gatk/
 * Picard Tools : https://broadinstitute.github.io/picard/
 * BBTools : https://sourceforge.net/projects/bbmap/
 * PrinSeq : http://prinseq.sourceforge.net/

## Authors
- [Charles Sanfiorenzo Cruz]() - charles.sanfiorenzo@upr.edu
- [Jenelys Ruiz Ortiz]() - jenelys.ruiz@upr.edu


## Contributors
 
  

## Disclaimer
This program takes advantage of a lot of software tools developed by many academic institutions. Make sure to cite their use in accordance to their policies!

## Citation
Sanfiorenzo C. J. & Ruiz J. (2016). MultiTool - A Bioinformatic Toolset [Computer Software]. University of Puerto Rico 
Rio Piedras. Retrieved from: https://github.com/CharlesSanfiorenzo/Bioinformatics.git

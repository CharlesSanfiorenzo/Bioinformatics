# Optimal ORF Finder (OOF)

## What does OOF do?
OOF takes in a fasta file containing one or more sequences and identifies all nested Open Reading Frames (ORFs) found within each sequence.

## Why was OOF created?

### Gene splicing
Genes often undergo gene splicing, a post-transcriptional modification that permits a single gene to code for multiple proteins. Gene splicing events occur through differential inclusion or exclusion of regions of pre-mRNA, producing protein isoforms that are often structurally and functionally distinct. But what about the stuff that gets spliced out? When splicing gets carried out in some genes, nested ORFs become available for translation - potentially producing short peptide products with unknown functions!

#### RNA cleavage
Much like splicing, RNA cleavage in the cytosol can lead to new, accessible ORFs. It is important to note, however, that the fragments possessing the new ORF(s) won't be capped. 

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Requirements](#requirements)
- [Support](#support)
 - [Help](#configuration-issueshelp)
 - [Bugs](#bugs--issues)
 - [Feature Requests](#feature-requests)
 - [Pull Requests](#pull-requests)
- [Authors](#authors) 

## Installation
OOF is entirely written in Python, meaning that all modules can be used once Python and all packages listed in the requirements are installed. 

## Usage

```python OOF.py -f <fasta> -o <outName> [optional: -t]```
- OOF outputs two multifasta files containing: (i) all nested ORFs per header in input fasta, and (ii) the longest nested ORF per header in input fasta. Use -t flag to translate ORFs. Yup, it's that simple.

## Requirements
* **Languages**
  * Python 2.7 or Python 3+ : https://www.python.org/
    * Packages : BioPython
 * **If on Windows**
   * Linux shell emulator : https://www.cygwin.com/
   * GNU Packages : http://gnuwin32.sourceforge.net/
  
## Support

### Configuration issues/help
If you need any help please contact one of the [authors](#authors) via email.

### [Bugs / Issues](https://github.com/CharlesSanfiorenzo/Bioinformatics/issues)
If you discover a bug in the script, please [create a new issue](https://github.com/CharlesSanfiorenzo/Bioinformatics/issues/new).

### [Feature Requests](https://github.com/CharlesSanfiorenzo/Bioinformatics/labels/Feature%20Request)
If you have a feature you would like added, please [create a new request](https://github.com/CharlesSanfiorenzo/Bioinformatics/issues/new).

### [Pull Requests]()
If you'd like to make your own changes ensure your PR is made against the ['dev']() branch.

## Authors
- [Charles Sanfiorenzo Cruz](https://github.com/CharlesSanfiorenzo/) - charles.sanfiorenzo@upr.edu | csanfior@mit.edu

## Citation
Sanfiorenzo C. (2017). Optimal ORF Finder (OOF) [Computer Software]. Massachusetts Institute of Technology. Retrieved from: https://github.com/CharlesSanfiorenzo/Bioinformatics.git

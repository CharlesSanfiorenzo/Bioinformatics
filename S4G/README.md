# Statistics for Genomes (S4G)

## What does S4G do?
S4G was a module originally part of ![MultiTool](https://github.com/CharlesSanfiorenzo/Bioinformatics/tree/master/MultiTool) that produced simple statistics for fully annotated genomes. Given the common usage of these metrics in genome-wide studies, S4G was repurposed into a fully functional, stand-allone shell script. Currently, S4G calculates **genome size**, **total intergenic length**, **average intron size**, **total intron length**, **count of intron-containing genes**, **intron count**, **CpG dinucleotide count**, and **GC Content**. 

    Total genes with introns: 55553
    Total introns: 62848
    Total intron length: 330523
    Average intron length: 5259.00902962
    Intergenic length: 91200324
    Genome size: 2906061972
    Intergenic/Genome Ratio: 0.0314
    CpG sites: 20052674
    GC Content: 24%

S4G also produces bed annotation files that contain intron positions per gene (which are not annotated in most genome annotation files found online).

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
Once Bedtools, Python, and all packages listed in the requirements are installed,  simply run  ```chmod +x S4G.sh``` to make it executable. Note that S4G creates & runs temporary python scripts in order to carry out some of its processes.

## Usage

```./S4G.sh <genome.fna> <genome.gbff> <genome.gff> > GenomeStats.txt```
- S4G will print out a series of statistics (see above) pertaining the genome annotations that were provided to it.

## Requirements
* **Languages**
  * Python 2.7 : https://www.python.org/
    * Packages : Pandas
* **If on Windows**
   * Linux shell emulator : https://www.cygwin.com/
   * GNU Packages : http://gnuwin32.sourceforge.net/
* **Tools**
   * Bedtools : http://bedtools.readthedocs.io/en/latest/
  
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
Sanfiorenzo C. (2017). Statistics for Genomes (S4G) [Computer Software]. University of Puerto Rico - Rio Piedras. Retrieved from: https://github.com/CharlesSanfiorenzo/Bioinformatics.git


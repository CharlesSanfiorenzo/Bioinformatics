# Motif Effect Predictor (MEP)

## Concept

However, in order to calculate **translational** efficiency one must rule out transcriptional regulation as a cause for reported change in translation. Commonly studies use the following formula to define translational regulation: . Because Ribosome profiling and RNA-Seq data sets are different in nature and yield different results when properly quantified, the aforementioned formula thus compares two incomparable values - it is comparing apples and oranges. To bypass this issue, MEP uses a simple equation that implements a retention value (r) to rule out transcriptional regulation in samples when quantifying translational efficiency:.

## Table of Contents
- [Installation](#installation)
- [Tools](#tools)
- [Requirements](#requirements)
- [Support](#support)
 - [Help](#configuration-issueshelp)
 - [Bugs](#bugs--issues)
 - [Feature Requests](#feature-requests)
 - [Pull Requests](#pull-requests)
- [Authors](#authors) 

## Installation
MEP is entirely written in Python, meaning that all modules can be used once Python and all packages listed in the requirements are installed. 

## Tools
>Note: AutoAligner.sh can be run to perform automated TopHat alignments, while AutoExpression.sh will perform automated >Cufflinks expression analysis and parse resulting files into the required format for MEP. In a scenario where >AutoExpression.sh is not able to find the Cufflinks installation path, one may perform the following GNU operations on the >resulting 'genes.fpkm_tracking' file produced by Cufflinks:
>
>```awk '{ print $1, $10, $14, $18, $22}' genes.fpkm_tracking | awk '{if (NR!=1) {print}}' > GenesFPKM.txt```

### Transcriptional Efficiency Estimator
Calculates transcriptional efficiency by comparing FPKM or RPKM values reported for genes. For transcriptional efficiency to be calculated, two RNA-Seq data sets are required in the Cufflinks expression analysis.

```python TranscriptionEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -R <regulation type (default:down)> ] > MotifsNoRef.fa```
 
### Translational Efficiency Estimator
Calculates translational efficiency by comparing FPKM or RPKM values reported for genes. For translational efficiency to be calculated, MEP requires that the user includes two Ribosome profiling data sets and two RNA-Seq data sets in the Cufflinks expression analysis. The retention value (r) is defined to be 0.30 by default.

```python TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type (default:down)> ] > MotifsNoRef.fa```
 
### MEP Database Creator
After calculating Transcriptional and/or Translational Efficiency, running this script will create a database for MEP's plotting & SVM modules. The resulting file (MEP.db) will be a csv containing the following info: Identifier, Motif Sequence, AlignmentS core, and Transcriptional/Translational Efficiency. More features can be added to MEP's data base, see documentation for more details.

```python DBCreator.py -f <fasta> -M <MEP.fpkm.csv> -N <NeedleResults>```

- #### Position Dependent Alignment (PDA)
  If the user is interested in the identity of bases in specific nucleotide positions across the motifs found in the multifasta file, a position dependent alignment score can be provided as an alternative to Needle scores in MEP's DBCreator as follows:
  
  ```python PDA.py -m <PWM Matrix> -f <fasta> > AlignmentScores.txt```
  
  ```python DBCreator.py -f <fasta> -M <MEP.fpkm.csv> -N AlignmentScores.txt```
  
  Note that the PWM has to be created by the user. Here is an example of a PWM for an 11 nucleotide motif:
  
                       0, 1, 0, 0
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 1, 0, 0
                       0, 1, 0, 0
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 1, 0, 0
                       0, 0, 1, 0
                       
  Column order represents: A, C, G, T. The PWM should be saved in a text file and provided to PDA.py.
  
### Relationship Plots
  MEP is capable of producing plots for the following scenarios:

- #### Multivariate Relationships
![](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MEP/doc/images/MultiVariate.png?raw=true)
- #### Single Feature Relationships
![](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MEP/doc/images/SingleFeature.png?raw=true)

### Effect Classifier


Please read the documentation for all available features and sub-commands, as well as in-depth guidance.

## Requirements
* **Languages**
  * Python 2.7 or Python 3+ : https://www.python.org/
    * Packages : Pandas, Seaborn, Scikit-learn
 * **If on Windows**
   * Linux shell emulator : https://www.cygwin.com/
   * GNU Packages : http://gnuwin32.sourceforge.net/
* **Tools**
  * Bowtie2 : http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  * TopHat : https://ccb.jhu.edu/software/tophat/manual.shtml
  * SAMtools : http://samtools.sourceforge.net/
  * Cufflinks : http://cole-trapnell-lab.github.io/cufflinks/
  
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
 
## Disclaimer
This program takes advantage of a lot of software tools developed by many academic institutions. Make sure to cite their use in accordance to their policies!

## Citation
Sanfiorenzo C. (2017). Motif Effect Predictor (MEP) [Computer Software]. Massachusetts Institute of Technology. Retrieved from: https://github.com/CharlesSanfiorenzo/Bioinformatics.git

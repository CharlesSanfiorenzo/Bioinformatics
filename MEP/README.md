# Motif Effect Predictor (MEP)

## Concept
DNA, RNA, and Protein sequences all possess embedded motifs that ascribe functions in cells. These functions can be regulatory, or they can be structural - both of which can, depending on the identity of the motif sequence, have different degrees of functionality. Motifs possess several hallmarks within or around the motif sequence itself, and are often thus classified in groups regarding their composition and function. MEP aims at elucidating the regulatory potential of a given group or class of motifs by integrating sequence information with results from the expression analysis of RNA-Seq and/or Ribosome Profiling data sets. If a motif is said to regulate gene expression at the level of **transcription**, for example, a user will want to have a measure of the efficiency of transcription for sequences that contain or are targetted by said motif, where <img width="350" height="41" src="https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MEP/doc/images/TranscEff.png?raw=true">. However, in order to calculate **translational** efficiency (for genes regulated at the level of translation) one must rule out transcriptional regulation as a cause for reported change in translation. Most studies use the following to define translational efficiency: ![](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MEP/doc/images/OldTE.png?raw=true). Because Ribosome profiling and RNA-Seq data sets are different in nature and yield different results when properly quantified, the aforementioned formula thus compares two incomparable values - it is comparing apples and oranges. To bypass this issue, MEP uses a simple equation that compares Ribosome profiling data sets by implementing a retention value (r) to rule out transcriptional regulation in samples when quantifying translational efficiency (read below). Transcriptional and/or Translational Efficiency values are then integrated into a multivariate analysis alongside sequence-related features (e.g. alignment score when comparing to a 'reference' or 'ideal' motif given by the user or inferred by MEP, unfolding energy, composition statistics, distance of motif to genic or genomic regions, etc) with the goal of comparing correlation values. Once a good relationship has been found, a quantified value of a given motif's regulatory effect can be predicted through simple linear regresion. This is particularly useful for engineering one's own motifs based on the degree of regulation one wishes to ascribe. If MEP shows correlation between multiple sequence-related features and Transcriptional or Translational Efficiency, two or more features can be used to infer if regulation occurs or not given a set of new, unclassified motifs throughout MEP's Support Vector Machine module.

**Use MEP if all of the following apply: (i) you are interested in quantifying transcriptional/translational efficiency of several genes that possess or are targetted by a group of motifs, (ii) possess RNA-Seq and Ribosome Profiling data sets in order to calculate transcriptional/translational efficiency, (iii) possess a multifasta file containing the motifs that affect the genes under study.**

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
>**Note: AutoAligner.sh can be run to perform automated TopHat alignments, while AutoExpression.sh will perform automated   Cufflinks expression analysis and parse resulting files into the required format for MEP. In a scenario where  AutoExpression.sh is not able to find the Cufflinks installation path, one may perform the following GNU operations on the  resulting 'genes.fpkm_tracking' file produced by Cufflinks:**
>
>```awk '{ print $1, $10, $14, $18, $22}' genes.fpkm_tracking | awk '{if (NR!=1) {print}}' > GenesFPKM.txt```

### Transcriptional Efficiency Estimator
Calculates transcriptional efficiency by comparing FPKM or RPKM values reported for genes. For transcriptional efficiency to be calculated, two RNA-Seq data sets are required in the Cufflinks expression analysis.

```python TranscriptionEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -R <regulation type (default:down)> ] > MotifsNoRef.fa```
 
### Translational Efficiency Estimator
Calculates translational efficiency by comparing FPKM or RPKM values reported for genes. For translational efficiency to be calculated, MEP requires that the user includes two Ribosome profiling data sets and two RNA-Seq data sets in the Cufflinks expression analysis. The following equation is used to calculate translational efficiency:

<p align="center">
  <img src="https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MEP/doc/images/TransEff.png?raw=true">
</p>

The retention value (r) is defined to be 0.30 by default. 
-R dictates the type of regulation that is being studied (i.e. up or down). The formulas above reflect downregulation. For up regulation studies, the efficiency ratio is calculated using the Experimental/Treated data sets as the denominator.

```python TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type (default:down)> ] > MotifsNoRef.fa```
 
### MEP Database Creator
After calculating Transcriptional and/or Translational Efficiency, running this script will create a database for MEP's plotting & SVM modules. The resulting file (MEP.db) will be a csv containing the following info: Identifier, Motif Sequence, Alignment Score, and Transcriptional/Translational Efficiency. More features can be added to MEP's data base, see documentation for more details.

```python DBCreator.py -f <fasta> -M <MEP.fpkm.csv> -N <NeedleResults>```

- #### Position Dependent Alignment (PDA)
  If the user is interested in the identity of bases in specific nucleotide positions across the motifs found in the multifasta file, a position dependent alignment score can be provided as an alternative to Needle scores in MEP's DBCreator as follows:
  
  ```python PDA.py -m <PWM Matrix> -f <fasta> > AlignmentScores.txt```
  
  ```python DBCreator.py -f <fasta> -M <MEP.fpkm.csv> -N AlignmentScores.txt```
  
  Note that the PWM has to be created by the user. Here is an example of a PWM for the following 11 nucleotide motif:
  **5'-CTTTCCTTTCG-3'**
  
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
These plots are produced from all the features stored in MEP's database files. Said features can be generated automatically by MEP (like Transcriptional/Translational Efficiency and Alignment Score), and others need to be added manually by the user (```DBCreator.py ``` has an option for this). I have been working on an automatic, wide-scale screening algorithm in order to autonomously provide a variety of metrics that have been proven in the literature to influence transcriptional/translational efficiency. Feel free to email [me](#authors) for updates.
- #### Single Feature Relationships
![](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MEP/doc/images/RiboBind.png?raw=true)

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

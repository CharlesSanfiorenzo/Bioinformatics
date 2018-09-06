# About
This repository contains several projects that have been funded by Puerto Rico's Louis Strokes Alliance for Minority Participation (PR-LSAMP) and MIT's Summer Research Program (MSRP). All projects within this repository are open-source, and may be used and modified as needed under educational purposes.


# Projects

## Motif Effect Predictor (MEP)
MEP is a set of computational tools capable of discovering and predicting relationships between several sequence-related features of DNA/RNA motifs, and transcriptional/translational efficiency. MEP currently evaluates features such as alignment score (when the motif under study is being compared to an 'ideal' motif), unfolding energy in areas surrounding the motif, statistics of the motif's composition of bases, and the motif's distance to other genomic or genic regions. All of the aforementioned metrics are compared to transcriptional and/or translational efficiency values derived from RNA-Seq and Ribosome Profiling data sets.
![](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MEP/doc/images/RiboBind.png?raw=true)

## MultiTool
MultiTool is an interactive multi-purpose bioinformatic tool that enables researchers with no background in bioinformatics to complete common tasks in the ever-growing field of '-omics'.

How it operates:

    $ ./MultiTool.sh 
    1) Extract         4) ReadQuality    7) BaseCounter   10) PSMC
    2) GeneCopies      5) ReadFilter     8) SNPCalling    11) Help
    3) FindGene        6) ProteomeSize   9) MutationLoad  12) Quit
    Enter your choice (use number):
The user can navigate through a series of curated options without the need of complex pipelines.

## Optimal ORF Finder (OOF)
How many nested Open Reading Frames (ORFs) can you find in a given genome, contig, chromosome, or gene? Run OOF to find out! Want to translate all of these nested ORFs and study their amino acid sequences? No problem, OOF has you covered! Need to identify the longest ORF within a set of nested ORFs? Say no more, OOF is the answer for you. 
<p align="center">
    <img width="550" alt="portfolio_view" src="https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/OOF/doc/OOF.png">
</p>

## Statistics for Genomes (S4G)
S4G is a hybrid shell/python module that calculates **genome size**, **total intergenic length**, **average intron size**, **total intron length**, **count of intron-containing genes**, **intron count**, **CpG dinucleotide count**, and **GC Content** of fully annotated genomes.
    
    $ ./S4G.sh genome.fa genome.gbff genome.gff
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

## Thermodynamic Promoter Activity Predictor (TPAP) 
TPAP is a DNA-protein binding prediction model for RNA Polymerase-promoter association that takes into consideration the thermodynamic principles behind TF-promoter binding events. It uses a supervised, binary Steiner tree model to predict changes in promoter activity when TF binding regions of a promoter are mutated or completely inhibited. 

## Main Author
- [Charles Sanfiorenzo](https://github.com/CharlesSanfiorenzo) - csanfior@caltech.edu | csanfior@mit.edu | charles.sanfiorenzo@upr.edu


## Contributors
- [Jenelys Ruiz Ortiz]() - jenelys.ruiz@upr.edu

## Disclaimer
Some of these projects take advantage of a lot of software tools developed by many academic institutions. Make sure to cite their use in accordance to their policies!

# Help Manual
Note: Try to run this script in a screen.
## I. Gene Extraction
This feature makes use of NCBI's Genbank annotated format files (GBKs) to produce a multifasta file containing all Open Reading Frames (ORFs) belonging to protein coding genes. First collect genbank files belonging to an organism's chromosomes (if assigned). Run MultiTool.sh w/ Extract option - this will output a file containing all ORFs for all protein coding genes. You can then choose to extract all ORFs or to only keep the longest ORF per gene (canonical sequences).

## II. Gene Copy Number Analysis
The majority of genes only possess a single functional copy within the genome of most organisms. We can take advantage of this fact to estimate the number of copies for a specific gene by calculating Whole-Genome Sequence (WGS) coverage for said gene relative to the average WGS coverage for all annotated genes across the genome. To do this, you would need to have sequence files containing the longest ORF per annotated gene/transcript (i.e. their canonical sequence, which contains all possible exons and introns), and fastq-formatted sequencing reads belonging to WGS projects. MultiTool can provide the user with canonical gene sequences given a genbank annotation file by using the Gene Extraction option (Note: choose 'Longest Transcript' when prompted).

### General Procedure

- **A.** Run MultiTool.sh and choose 'Extract' option. Choose 'Longest Transcript' when prompted. 
    * **Output:** All canonical ORFs in fasta format.
- **B.** Run MultiTool.sh and choose 'Find Gene' option. Input the gene symbol of the gene that's going to be quantified. 
    * **Output:** Gene under study in fasta format. 
- **C.** Run MultiTool.sh and choose 'Copy Number' option. If the user needs to quantify multiple genes at once, or has yet to quantify the average coverage for all canonical ORFs, MultiTool.sh has an option for this (choose 'MultipleAlignments' when prompted, else just choose 'Single Alignment'). The reference files used should contain the gene sequence under study in fasta format (produced in Step B), and all canonical ORFs in fasta format (produced in Step A). 
    * **Output:** Copy Number of gene(s) selected.

### What is MultiTool doing in order to assess gene copy number?
MultiTool is simply using Bowtie2 & SAMtools in order to perform a series of pre-processing, alignment, and post-processing operations. The most important aspects of MultiTool's Gene Copy Number functionality are the alignment and post-processing steps:
- Alignment: MultiTool uses a modified, stringent version of Bowtie2's standard commandline.

```bowtie2 -x $Build -D 0 -R 0 -N 0 --no-1mm-upfront --no-mixed --no-discordant $Reads --threads 16 | samtools view -bS - > $Build.bam```

For every read segment, no mismatches or mixed-segment alignments are allowed. Segments that only align partially to a region are also discarded, and no re-seeding is performed. This allows MultiTool to only quantify reads that fully align to a particular region.
- Post-Processing: To take care of over-sampling or "over-representation" of the genome, MultiTool uses PicardTools to target duplicate alignments (i.e. alignments that correspond to multiple reads mapping to the same locations) and only keep the one with the highest score.

```java -XX:ParallelGCThreads=16 -jar picard.jar MarkDuplicates INPUT=BuildName.sorted.bam OUTPUT=NoDup.bam METRICS_FILE=metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT```

### Pitfalls
While this approach works for most genes, genes that have undergone retro-insertion events in the genome (i.e. genes that possess retro-copies) tend to cause MultiTool to integrate many false positives into its analysis. Retro-copies have an increasingly high number of frameshift mutations, potentially making the product of the retro-gene in question non-functional in respect to the cellular functionality of its original copy.

## III. Find Gene
This feature extracts a specific gene header and sequence from a mutifasta file. Run this script w/ Extract option. When prompted, use the fasta output from the Extract option as input; specify the gene of interest. Inclusive mode will keep all genes with similar names, and exclusive mode will keep the gene named with the regular expression given.

## IV. Read Quality Test for WGS
Using fastq files & a reference fasta, this option will output alignment rate, average coverage, and Transition/Transversion ratio (Ts/Tv) for the selected organism. 

Note: Only FastQC screening is needed to test for quality for tNGS methods, as multiple sequencing runs over the same target increase validity of base calls.

### Walkthrough

Reads are the raw output of a sequencing machine. They need to be tested for quality, in order to determine whether re-sequencing is necessary. Before using this option, reads should be screened with the use of FastQC. If evaluating WGS, GC content of the reads should be close or equal to GC content in the genome assembly and should follow a Gaussian distribution. All reads should follow Chargaff's rule (seen in content per base section of FastQC), as strands have a 1:1 sequencing probability. Average Quality should always be above 20, and length per read should follow a Triangular distribution unless they were trimmed post-sequencing. If sequenced reads have a decent coverage fold of the genome, filtering reads is preferred over trimming.

[This](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html#M0) is how good quality reads look like in general when screened with FastQC.

However, for WGS, screening with FastQC is not enough. Once we are sure we have what appear to be good quality reads, we must perform a second screening after alignment steps...

Good Quality indicators for WGS mammalian data (post-alignment)

Attribute | Results
------------ | -------------
Alignment Rate | More than 90%
Average Coverage| More than 12
Ts/Tv Ratio | Higher than or close to 2.1

For all whole-genome studies, no matter the class of the sampled organism, Ts/Tv should be close to or above 2.0. Values below 1 seem to be related to the introduction of random sequencing errors, which bypass initial FastQC screenings. These errors can be observed when binning the Ts/Tv over the whole genome. Also, random sequencing errors greatly increase the number of indels reported, and reduce the number of SNPs reported (there is a reason for this, see Part VII, section G below). Protein coding percentage and CpG Island spans in the genome positively correlate with Ts/Tv (this is due to the often methylated cytosines in CpG dinucleotides and overall probability of isochore formation), so there should be a trend when comparing the Ts/Tv's of different species.

Another indicator of good quality is the pKaKs ratio. Random sequencing errors skew the value for pKaKs close to 1 (i.e. neutral), which is not possible. A pKaKs ratio of 1 implies that no change is occurring at a species level (i.e. it goes against the concepts behind evolution). All species have a fixed pKaKs, even when selective pressures are added into the mix. Hence, assessing and comparing this value for the organism that one is studying provides an additional quality assessment layer for WGS. 

## V. Proteome Size
This feature will annotate the amount of aminoacids or codons in an organism's proteome by using sequences stored in Genbank or RefSeq files. It will also count the headers (amount of protein coding genes annotated).

1. Download from NCBI's ftp page a genbank belonging to either all RNA or protein sequences.
2. Choose which type of data you are using (i.e. transcriptomic vs proteomic)
3. Verify output: 'ProteomeSize.txt'.

## VI. Base Counter
Outputs amount of bases in individual or multifasta files. If you wish to determine the amount of bases per sequence in a multifasta file, please use 'counter.pl'.
* Use: perl counter.pl inputfile.fa > outputfile

## VII. SNP Calling
Using fastq files & a reference fasta genome, this feature produces a quality-filtered VCF and a statistics file containing information related to the amount of SNPs and Indels in the organisms genome. You may choose to run Bowtie2 and Samtools mpileup alone, or cross-analyze your data using Bowtie2/GATK, Bwa-mem/GATK, Bowtie2/Samtools, and Bwa-mem/Samtools combinations.
* Input: fastq files, reference genome assembly (fasta)
* Output: Filtered VCF, VCF statistics file

### Walkthrough 

If one wishes to determine single nucleotide polymorphisms (SNPs) within a particular region of a genome (e.g. a gene), or accross the whole genome itself, having a curated SNP Calling pipeline is ultimately necessary for producing accurate inferences. However, the amount and identity of the SNPs called depend heavily on the quality of your data set. If not sure how to test read quality, please refer to part IV. The following snippet will illustrate how a curated SNP calling pipeline looks, in accordance to [Broad Institute's guidelines](https://software.broadinstitute.org/gatk/best-practices/) for WGS :

    # A. Index genome to a Bowtie2 build (provides fast access to genome information)
    bowtie2-build GenomeAssembly.fna BuildName
    # B. Align sequenced data (reads) to build using Bowtie2
    bowtie2 -x BuildName -1 read_1 -2 read_2 --threads 16 -S BuildName.sam
    # C. Convert from sam format to bam format (bam is the binary version of a bam -> occupies less space, fast random access)
    samtools view -bS BuildName.sam -@ 16 > BuildName.bam
    # D. Sort bam file concordantly using SAMtools
    samtools sort BuildName.bam BuildName.sorted -@ 16
    # E. Mark and remove duplicates in the alignment using Picardtools
    java -XX:ParallelGCThreads=16 -jar picard.jar MarkDuplicates INPUT=BuildName.sorted.bam OUTPUT=NoDup.bam METRICS_FILE=metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    # F. Define indel intervals using GATK (for realignment, seen in next step)
    java -XX:ParallelGCThreads=16 -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R Genome.fna -I NoDup.bam -o Intervals.txt
    # G. Realign against indels using GATK
    java -XX:ParallelGCThreads=16 -jar GenomeAnalysisTK.jar -T IndelRealigner -R Genome.fa -I NoDup.bam -targetIntervals Intervals.txt -o NoDupRealigned.bam
    # H. Call SNPs (variants) using SAMtools mpileup
    samtools mpileup -uf Genome.fa NoDupRealigned.bam | bcftools call -c > BuildName.vcf
    # I. Filter VCF based on BAQ or SNP depth using VCFtools
    vcftools --vcf BuildName.vcf --recode --recode-INFO-all --out BuildName --minQ 20.00 
    # J. Produce statistics for the VCF 
    bcftools stats BuildName.recode.vcf > BuildNameFilteredVCFStats.txt
Before explaining each of the steps listed in the code above, it's important to understand the software options available for alignments and variant calling. We are using Bowtie2 as our aligner and Samtools mpileup as our variant caller in this example, but a large group of researchers tend to use Bwa-mem for alignments and GATK's Haplotype caller to call for variants. The decision to choose between Bowtie2 and Bwa-mem, or Samtools mpileup and GATK's Haplotype caller is arbitrary, although GATK standarized variant calls are optomized for human genomes. Nonetheless, Bowtie2 and Samtools mpileup mash-ups have been favored more often in recent publications. If unsure, this script allows one to run alignments and variant calls combinatorially for comparison purposes.

#### A. Index genome
Most aligners require indexing of the reference sequence in order to increase aligning and processing speed. The reference sequence belongs to a healthy individual, and is used as the norm for sequence comparison. This is a standard procedure for Bowtie2 and Bwa-mem.
#### B. Alignment
Alignment between sequenced data and a reference sequence. Usually a genome assembly is used as the reference sequence, as chromosome scaffolds or individual gene sequences lead to biased variant calls due to forced alignments. Default parameters in most alignment tools are ideal, unless the intended analysis consists of comparing sequences belonging to two different species; localized and sensitive parameters are preferred for the latter.
#### C. Sam to bam conversion
This step can be integrated into the alignment if needed. Having a binary and compressed version of the alignment improves processing speed and increases space availability. 
#### D. Sorting a bam file concordantly
Most downstream analysis tools require for bam files to be sorted. This allows easier access to compressed annotations.
#### E. Marking and removing duplicates
In WGS, sometimes reads map multiple times to a single region in the genome. This is particularly true for sequences that were amplified through PCR using low complexity libraries. We can remove duplicate alignments with Picardtools; average coverage fold of reads should be re-calculated again after this step if evaluating a whole genome.

Note: This step should be skipped for RNAseq studies, as it is expected for reads to map multiple times to particular reference sequences (due to expression levels). 

#### F. Determine indel intervals for realignment
All algorithms used for alignment have problems with aligned sequences close to regions identified as possible indels. To improve the accuracy of our calls in Step H, coordinates for indel sites are annotated in this step and reassessed in the next step.

Note: Indel realignment steps may be skipped when using GATK's Haplotype caller or any variant caller that performs a haplotype assembly step.

#### G. Realignment around the indels
Using the indel coordinates produced in the previous step, GATK realigns sequences against indels to determine if the annotated indel is truly an indel, in this way fixing alignments that may have been originally considered mismatches accross large regions of the reads. If the data set belongs to an organism that has lists of known SNPs available, one may also perform an additional step consisting of base quality score recalibration (i.e. if the SNP is reported both in the list and the alignment, the certainty of the reported SNP increases).

#### H. Variant calling
SNPs and indels are annotated formally in this step. If a particular region of the genome was amplified, one may offer a set of coordinates for the variant caller to use. This lowers processing time exponentially.

Example:
     
    samtools mpileup -r 1:860-1000 -uf Genome.fa NoDupRealigned.bam | bcftools call -c > BuildName.vcf
    
where -r is the region consisting of chromosome#:bpLowerRange-bpUpperRange,
-u is uncompressed,
and -f is the reference file (e.g. assembly, canonical gene sequence, etc).

#### I. VCF Filtering
In order to remove calls with low certainty (low BAQ, low Q, and/or low SNP depth) from the analysis, we need to filter the VCF. Usually Q and BAQ values above 20 (less than 1/10<sup>2</sup> chance of being an incorrect call/alignment respectively) are excellent, so this value may be used as a minimum treshold. Some articles in the literature suggest that SNPs with depth over 1/3 of the average SNP depth be kept, but this is completely dependant on the type of sequencing and analysis one is doing.

Oh, and for future reference:

- Q = Phred quality score (certainty that a base call during sequencing is correct)
- BAQ = Base alignment quality (certainty that an alignment between two bp is correct, depends on neighboring sequences)

#### Bonus Step
One can assess the possible phenotypic impact of SNPs in coding and regulatory regions by running snpEff using a VCF file. snpEff takes into account synonymous and non-synonymous changes in codons, as well as similarity and dissimilarity between the resulting amino acid residues if the gene is protein coding. Loss or addition of start and stop codons is also evaluated.

#### One final note: 
Targeted Next-Generation Sequencing (tNGS) is considered a fast & cost-effective method to detect multiple mutations simultaneously in regions of interest. While it may seem like a good idea for studies involving particular mutations in regulatory regions or genes, the accuracy of the sequencing procedure may contribute to problems during alignment steps if using a genome assembly as reference rather than the canonical sequence of the targeted gene or regulatory region. The more one knows about the sequences around our target, the better we can identify it visually - the same applies for alignment algorithms. A caveat around this would be cloning our target sequence, doing a whole-genome sequencing run with our sample, doing multiple tNGS runs with our target, and then comparing SNP statistics within that region and analyzing the alignments visually using IGV.

## VIII. Mutational Load
This feature assays the amount of possible synonymous and non-synonymous sites in a gene's canonical coding sequence, as well as the amount of synonymous and non-synonymous SNPs in the coding regions of raw sequencing runs (we use snpEff for this last bit). These values are used to calculate the polymorphic Ka/Ks ratio of the genome through the Nei-Gojobori method, which is described as follows:
![equation](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MultiTool/docs/images/NeiGojobori.png?raw=true)

* Input: Filtered VCF, snpEff Library
* Output: SNPStats.txt (SNP site and occurance info), pKa/Ks.txt (polymorphic Ka/Ks)

The output also includes a theoretical value of the average impact of SNPs across protein-coding regions of the genome derived from classifications in snpEff. The average is calculated in the following way:
![equation](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MultiTool/docs/images/MeanSNPImpact.png?raw=true)

## IX. PSMC
[PSMC](https://github.com/lh3/psmc) is a software package capable of inferring population size history from diploid genomic sequences using the Pairwise Sequentially Markovian Coalescent (PSMC) model. This feature provides a working pipeline capable of running PSMC from multiple vcf files from different species. However, amongst the parameters needed for a correct population size inference, mutation rate(u), generation times (g), and a factor for the evolutionary distance to the most recent common ancestor (t) need to be changed for each species under study. If u, g, and t are different for a,b,c..z species, PSMC's default operations are not capable of scaling and comparing population histories for said species. This package gives users the option to input multiple n's, g's,and t's in accordance to the organisms being studied and compared - resulting in a plot containing traces of the population history pertinent to each species. In this example, a single mutation rate (u) is used for simulated data representing species with multiple generation times (g's, found next to species name).
![](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MultiTool/docs/images/psmcExample.png?raw=true)

Note that the values for u and g can be found in the literature. However, while the model is capable of estimating a value for t per loci, the theoretical maximum value for t must be inferred by the user. MultiTool includes tools capable of analyzing the output of PSMC by showing the distribution of t's estimated by the model. 

![](https://github.com/CharlesSanfiorenzo/Bioinformatics/blob/master/MultiTool/docs/images/FigureTs.png?raw=true)

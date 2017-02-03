#Help Manual
Note: Try to run this script in a screen.
##I. Gene Extraction
This feature makes use of NCBI's Genbank annotated format files (GBKs) to produce a multifasta file containing all ORF's belonging to protein coding genes. First collect genbank files belonging to an organism's chromosomes (if assigned). Run MultiTool.sh w/ Extract option; this will output a file containing all ORF's for all protein coding genes. You can choose to extract all ORFs or to only keep the longest ORF per gene (canonical sequences).

##II. Gene Copy Number Analysis
Run this script w/ Extract option, but keep only canonical sequences. Parse the output to extract desired genes using 'Find Gene'. Then, map reads to individual genes, and then to all ORF's using GeneCopies option. You should end up with a coverage result for the individual gene in question (Ref: GeneName.fa), and for all of the ORFs. Divide the coverage value belonging to each individual gene between the coverage value for all ORFs (i.e. gene/AllORFs). Since most genes have a single copy, this result should be an estimate value for the copy number of that gene.

##III. Find Gene
This feature extracts a specific gene header and sequence from a mutifasta file. Run this script w/ Extract option. When prompted, use the fasta output from the Extract option as input; specify the gene of interest. Inclusive mode will keep all genes with similar names, and exclusive mode will keep the gene named with the regular expression given.

##IV. Read Quality Test
Using fastq files & a reference fasta, this option will output alignment rate, average coverage, and TS/TV ratio for the selected organism.

Good Quality indicators for genomic mammalian data

Atribute | Results
------------ | -------------
Alignment Rate | More than 90%
Average Coverage| More than 12
TS/TV Ratio | Higher than or close to 2.1


##V. Proteome Size
This feature will annotate the amount of aminoacids or codons in an organism's proteome by using sequences stored in Genbank or RefSeq files. It will also count the headers (amount of protein coding genes annotated).

1. Download from NCBI's ftp page a genbank belonging to either all RNA or protein sequences.
2. Choose which type of data you are using (i.e. transcriptomic vs proteomic)
3. Verify output: 'ProteomeSize.txt'.

##VI. Base Counter
Outputs amount of bases in individual or multifasta files. If you wish to determine the amount of bases per sequence in a multifasta file, please use 'counter.pl'.
* Use: perl counter.pl inputfile.fa > outputfile

##VII. SNP Calling
Using fastq files & a reference fasta genome, this feature produces a quality-filtered VCF and a statistics file containing information related to the amount of SNPs and Indels in the organisms genome. You may choose to run Bowtie2 and Samtools mpileup alone, or cross-analyze your data using Bowtie2/GATK, Bwa-mem/GATK, Bowtie2/Samtools, and Bwa-mem/Samtools combinations.
* Input: fastq files, reference genome assembly (fasta)
* Output: Filtered VCF, VCF statistics file

###Walkthrough 

If one wishes to determine single nucleotide polymorphisms (SNPs) within a particular region of a genome (e.g. a gene), or accross the whole genome itself, having a curated SNP Calling pipeline is ultimately necessary for producing accurate inferences. However, the amount and identity of the SNPs called depend heavily on the quality of your data set. If not sure how to test read quality, please refer to part IV. The following snippet will illustrate how a curated SNP calling pipeline looks, in accordance to Broad Institute's guidelines<sup>[1]</sup> for WGS :

    # A. Index genome to a Bowtie2 build (provides fast access to genome information)
    bowtie2-build GenomeAssembly.fna BuildName
    # B. Align sequenced data (reads) to build using Bowtie2
    bowtie2 -x BuildName -1 read_1 -2 read_2 --threads 16 -S BuildName.sam
    # C. Convert from sam format to bam format (bam is the binary version of a bam -> occupies less space, fast random access)
    samtools view -bS BuildName.sam -@ 16 > BuildName.bam
    # D. Sort bam file concordantly using SAMtools
    samtools sort BuildName.bam BuildName.sorted -@ 16
    # E. Mark and remove duplicates in the alignment using Piccardtools
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

####A. Index genome
Most aligners require indexing of the reference sequence in order to increase aligning and processing speed. Standard procedure for Bowtie2 and Bwa-mem.
####B. Alignment
Alignment between sequenced data and a reference sequence. Usually a genome assembly is used as the reference sequence, as chromosome scaffolds or individual gene sequences lead to biased variant calls due to forced alignments. Default parameters in most alignment tools are ideal, unless the intended analysis consists of comparing sequences belonging to two different species; localized and sensitive parameters are preferred for the latter.
####C. Sam to bam conversion
This step can be integrated into the alignment if needed. Having a binary and compressed version of the alignment improves processing speed and increases space availability. 
####D. Sorting a bam file concordantly
Most downstream analysis tools require for bam files to be sorted. This allows easier access to compressed annotations.
####E. Marking and removing duplicates


##VIII. Mutational Load
This feature assays the amount of possible synonymous and non-synonymous sites in a gene's canonical coding sequence, as well as the amount of synonymous and non-synonymous SNPs in the coding regions of the raw sequence of the organism (we use snpEff for this last bit). These values are used to calculate the polymorphic Ka/Ks ratio of the genome through the Nei-Gojobori method.
* Input: Filtered VCF, snpEff Library
* Output: SNPStats.txt (SNP site and occurance info), pKa/Ks.txt (polymorphic Ka/Ks)


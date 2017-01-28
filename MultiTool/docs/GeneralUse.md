#Help Manual
Note: Try to run this script in a screen.
#I. Gene Extraction
This feature makes use of NCBI's Genbank annotated format files (GBKs) to produce a multifasta file containing all ORF's belonging to protein coding genes. First collect genbank files belonging to an organism's chromosomes (if assigned). Run MultiTool.sh w/ Extract option; this will output a file containing all ORF's for all protein coding genes. You can choose to extract all ORFs or to only keep the longest ORF per gene (canonical sequences).

#II. Gene Copy Number Analysis
Run this script w/ Extract option, but keep only canonical sequences. Parse the output to extract desired genes using 'Find Gene'. Then, map reads to individual genes, and then to all ORF's using GeneCopies option. You should end up with a coverage result for the individual gene in question (Ref: GeneName.fa), and for all of the ORFs. Divide the coverage value belonging to each individual gene between the coverage value for all ORFs (i.e. gene/AllORFs). Since most genes have a single copy, this result should be an estimate value for the copy number of that gene.

#III. Find Gene
This feature extracts a specific gene header and sequence from a mutifasta file. Run this script w/ Extract option. When prompted, use the fasta output from the Extract option as input; specify the gene of interest. Inclusive mode will keep all genes with similar names, and exclusive mode will keep the gene named with the regular expression given.

#IV. Read Quality Test
Using fastq files & a reference fasta, this option will output alignment rate, average coverage, and TS/TV ratio for the selected organism.

Good Quality indicators for genomic mammalian data

Atribute | Results
------------ | -------------
Alignment Rate | More than 90%
Average Coverage| More than 12
TS/TV Ratio | Higher than or close to 2.1


#V. Proteome Size
This feature will annotate the amount of aminoacids or codons in an organism's proteome by using sequences stored in Genbank or RefSeq files. It will also count the headers (amount of protein coding genes annotated).

1. Download from NCBI's ftp page a genbank belonging to either all RNA or protein sequences.
2. Choose which type of data you are using (i.e. transcriptomic vs proteomic)
3. Verify output: 'ProteomeSize.txt'.

#VI. Base Counter
Outputs amount of bases in individual or multifasta files. If you wish to determine the amount of bases per sequence in a multifasta file, please use 'counter.pl'.
* Use: perl counter.pl inputfile.fa > outputfile

#VII. SNP Calling
Using fastq files & a reference fasta genome, this feature produces a quality-filtered VCF and a statistics file containing information related to the amount of SNPs and Indels in the organisms genome.
* Input: fastq files, reference genome assembly (fasta)
* Output: Filtered VCF, VCF statistics file

#VIII. Mutational Load
This feature assays the amount of possible synonymous and non-synonymous sites in a gene's canonical coding sequence, as well as the amount of synonymous and non-synonymous SNPs in the coding regions of the raw sequence of the organism (we use snpEff for this last bit). These values are used to calculate the polymorphic Ka/Ks ratio of the genome through the Nei-Gojobori method.
* Input: Filtered VCF, snpEff Library
* Output: SNPStats.txt (SNP site and occurance info), pKa/Ks.txt (polymorphic Ka/Ks)


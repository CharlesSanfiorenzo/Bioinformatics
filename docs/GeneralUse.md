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
2. Verify output: 'ProteomeSize.txt'.

#VI. Base Counter
Outputs amount of bases in individual or multifasta files. If you wish to determine the amount of bases per sequence in a multifasta file, please use 'Counter.pl'.
* Use:
perl Counter.pl inputfile.fa > outputfile

#VI. SNP Calling

         VII. PSMC Analysis

            File Naming Protocol (Examples):
               Reference: TP53.fa
               BuildName: TP53
               ResultsName: TP53Cov.txt
               etc.
               Other outputs will be named in accordance to the build.

        VIII. Tool Box
              A collection of tools for frequent tasks. Meant for
              bigginers in the terminal!
          

                 More will be added soon!
            
         Authors 
            Charles J. Sanfiorenzo Cruz 
            Jenelys Ruiz Ortiz

         Collaborators

         Version 1.4

#Help Manual
Note: Try to run this script in a screen.
#I. Gene Extraction
This feature makes use of NCBI's Genbank annotated format files (GBKs) to produce a multifasta file containing all ORF's belonging to protein coding genes. 
-1) Collect genbank files belonging to an organism's chromosomes (if assigned); if unassigned, use the single GBK available (Un). Concatenate all files into a single file.    
-2) Run MultiTool.sh w/ Extract option. This will output 'LongestGenome.fa', containing all ORF's for all protein coding genes.  

          II. Gene Copy Number Analysis
             -1) Run this script w/ Extract option.
             -2) Extract genes of interest from output (LongestGenome.fa).
             -3) Parse 'LongestGenome.fa' to extract desired genes using 
                 Find Gene.
             -4) Map reads to individual genes, and then to all ORF's using 
                 GeneCopies option. You should end up with a coverage result 
                 for the individual gene in question (Ref: GeneName.fa), and 
                 for all of the ORFs (Ref: LongestGenome.fa).
             -5) Divide the coverage value belonging to each individual gene 
                 between the coverage value for all ORFs (gene/AllORFs). This
                 is an estimate value for the copy number of that gene.

         III. Find Gene
              This feature extracts a specific gene header and sequence from a
              mutifasta file.
             -1) Run this script w/ Extract option.
             -2) When prompted, use the fasta output from the Extract option as
                 input; specify the gene of interest.
              Options:
               -Exclusive: Stringent; only extracts genes with name specified.
                Example
                    Sample= TP53_1, TP53_2, TP53_3, TP53_4, TP53RT2_1, TP52_1
                    Gene specified= TP53
                    Results= TP53_1
               -Inclusive: Extracts genes with similar names.
                Example
                    Sample= TP53_1, TP53_2, TP53_3, TP53_4, TP53RT2_1, TP52_1
                    Gene specified= TP53
                    Results= TP53_1, TP53_2, TP53_3, TP53_4, TP53RT2_1, TP52_1

         IV. Read Quality Test
             Using fastq files & reference fasta, this option will output
             alignment rate, average coverage, and TS/TV ratio for select
             organism.
             ->  Good Quality indicators for mammalian data: <-
              a) Alignment rate > 90%
              b) Average Coverage > 12
              c) TS/TV ratio higher than, or close to 2.1
            
          V. Proteome Size
               This feature will annotate the amount of aminoacids or codons in 
               an organism's proteome by using sequences stored in Genbank
               files. It will also count the headers (amount of protein coding
               genes annotated).
              -1) Download from NCBI's ftp page a genbank belonging to either 
                  all RNA or protein sequences. 
              -2) Verify output: 'ProteomeSize.txt'.

          VI. Base Counter
               Outputs amount of bases in individual or multifasta files.

           VI. SNP Calling

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

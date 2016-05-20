#!/bin/bash

#Dependencies: Bowtie2, SAMtools, BBTools, gb2tab, standard unix, python-2.7, perl, 

PS3='Enter your choice (use number): '
options=("Extract" "GeneCopies" "FindGene" "ReadQuality" "ProteomeSize" "SNPCalling" "Help" "Quit")
suboptions=("SingleRead" "ReadPair")
suboptionspsize=("RnaGBK" "ProteinGBK")

select opt in "${options[@]}"
do
    case $opt in
        "Extract")
echo "Type in the name of the GBK file"
read -p "GBK file : " gbk
echo "Processing..."
if [ -r $gbk ]; then
python ../../../../opt/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

cat << EOF > pyscript.py
#!/usr/bin/python

import fileinput

for line in fileinput.input():
	tokens = line.split("\t")
	name = tokens[0]
	seq = tokens[1]

	print ">%s\n%s" % (name, seq)

EOF

chmod +x pyscript.py

./pyscript.py Genomebaby.txt > Genome.txt

rm pyscript.py

. ../../../../opt/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. ../../../../opt/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
sed 's,__,%,g' Genomev3.fa | sed 's,_,+,g' | sed 's,+,_1,g' | sed 's,%,_,g' > TiddyGenome.fa
cat TiddyGenome.fa |\
awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' |\
tr "_" "\t" |\
awk -F '   '  '{printf("%s\t%d\n",$0,length($3));}' |\
sort -t '	' -k1,1 -k4,4nr |\
sort -t '	' -k1,1 -u -s |\
sed 's/	/_/' |\
cut -f 1,2 |\
tr "\t" "\n"  |\
fold -w 80 > LongestGenome.fa
echo "Process completed; please verify the outputs. Output can be found in 
      'LongestGenome.fa'."
rm Genome.txt Genome.fa Genomev2.fa Genomev3.fa TiddyGenome.fa Genomebaby.txt
exit 0
else 
echo "Input is missing or is in an incorrect format."
exit 1
fi
            ;;
        "GeneCopies")
clear
select subopt in "${suboptions[@]}"
do
    case $subopt in
        "SingleRead")
echo "Select your reference (.fa) and read (.fastq) files; include path."
read -p "Reference : " ref
read -p "Read : " reads
clear
echo "Select a name for your build and cov results."
read -p "BuildName :" build
read -p "ResultsName : " CovResults
cmd="sort"
echo $ref $build | bowtie2-build $ref $build
if [ -r $ref ]; then
echo $build $reads | bowtie2 -x $build -D 0 -R 0 -N 0 --no-mixed --no-discordant $reads --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo "Process complete. Please verify '$CovResults'."
exit 0
else 
echo "Outputs missing and/or truncated, please verify inputs"
exit 1
fi       
            ;;
        "ReadPair")
echo "Select your reference (.fa) and read (.fastq) files; include path."
read -p "Reference : " ref
read -p "Read_1 : " read_1
read -p "Read_2 : " read_2
clear
echo "Select a name for your build, sortedbam, and cov results."
read -p "BuildName :" build
read -p "ResultsName : " CovResults
cmd="sort"
echo $ref $build | bowtie2-build $ref $build
if [ -r $ref ]; then
echo $build $read_1 $read_2 | bowtie2 -x $build -D 0 -R 0 -N 0 --no-mixed --no-discordant -1 $read_1 -2 $read_2 --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo "Process complete. Please verify '$CovResults'."
exit 0
else 
echo "Outputs missing and/or truncated, please verify inputs"
exit 1
fi
;;
              *) echo invalid option;;
        
esac
done
           ;;
    "FindGene")
echo "Select your multifasta and desired gene to be extracted"
read -p "FastaFile : " multifasta
read -p "GeneName : " genename
clear
echo "Processing..."
if [ -r $multifasta ]; then
echo $multifasta $genename | awk 'BEGIN {RS=">"} /'$genename'/ {print ">"$0}' $multifasta | awk -v RS='>' 'NR>1 {gsub("\n", ";", $0); sub(";$", "", $0); print ">"$0}' | head -n 1 | tr ';' '\n' | sed '/^$/d' > $genename.fa
echo "Process complete. Output can be found in '$genename.fa'."
exit 0
else
echo "Input is missing and/or truncated"
exit 1
fi
           ;;
 "ReadQuality")
clear
select subopt in "${suboptions[@]}"
do
    case $subopt in
        "SingleRead")
echo "Select your reference (.fa) and read (.fastq) files; include path."
read -p "Reference : " ref
read -p "Read : " reads
clear
echo "Select a name for your build, vcf, cov & vcf results."
read -p "BuildName : " build
read -p "VcfName :" vcf
read -p "CovResultsName : " CovResults
read -p "VcfStatsName : " VcfResults
cmd="sort"
echo $ref $build | bowtie2-build $ref $build
if [ -r $ref ]; then
echo $build $reads | bowtie2 -x $build $reads --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $vcf
echo $vcf $VcfResults | bcftools stats $vcf > $VcfResults
exit 0
else 
echo "Outputs missing and/or truncated, please verify inputs"
exit 1
fi       
            ;;
        "ReadPair")
echo "Select your reference (.fa) and read (.fastq) files; include path."
read -p "Reference : " ref
read -p "Read_1 : " read_1
read -p "Read_2 : " read_2
clear
echo "Select a name for your build, sortedbam, vcf, cov & vcf results."
read -p "BuildName :" build
read -p "VcfName :" vcf
read -p "CovResultsName : " CovResults
read -p "VcfStatsName : " VcfResults
cmd="sort"
echo $ref $build | bowtie2-build $ref $build
if [ -r $ref ]; then
echo $build $read_1 $read_2 | bowtie2 -x $build -1 $read_1 -2 $read_2 --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $vcf
echo $vcf $VcfResults | bcftools stats $vcf > $VcfResults
exit 0
else 
echo "Outputs missing and/or truncated, please verify inputs"
exit 1
fi
;;
              *) echo invalid option;;
        
esac
done
           ;;
 "ProteomeSize")
clear
select subopt in "${suboptionspsize[@]}"
do
    case $subopt in
        "RnaGBK")
echo "Type in the name of the GBK file"
read -p "GBK file : " gbk
if [ -r $gbk ]; then
echo "Processing..."
python ../../../../opt/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

cat << EOF > pyscript.py
#!/usr/bin/python

import fileinput

for line in fileinput.input():
	tokens = line.split("\t")
	name = tokens[0]
	seq = tokens[1]

	print ">%s\n%s" % (name, seq)

EOF

chmod +x pyscript.py

./pyscript.py Genomebaby.txt > Genome.txt

rm pyscript.py



. ../../../../opt/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. ../../../../opt/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
sed 's,__,%,g' Genomev3.fa | sed 's,_,+,g' | sed 's,+,_1,g' | sed 's,%,_,g' > TiddyGenome.fa
cat TiddyGenome.fa |\
awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' |\
tr "_" "\t" |\
awk -F '   '  '{printf("%s\t%d\n",$0,length($3));}' |\
sort -t '	' -k1,1 -k4,4nr |\
sort -t '	' -k1,1 -u -s |\
sed 's/	/_/' |\
cut -f 1,2 |\
tr "\t" "\n"  |\
fold -w 80 > LongestProtein.fa
grep '^>' LongestProtein.fa | wc -l | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Genes"}' > GeneCount.txt
grep -v '^>' LongestProtein.fa | wc -c | awk '{print $1/3}' |  awk '{printf "%2.1f\n",$0}' | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Aminoacids"}' > AminoacidCount.txt
cat GeneCount.txt AminoacidCount.txt > ProteomeSize.txt
echo "Process completed; please verify the outputs."
echo "Proteome size can be found in 'ProteomeSize.txt'"
rm Genome.txt Genome.fa Genomev2.fa Genomev3.fa TiddyGenome.fa GeneCount.txt AminoacidCount.txt Genomebaby.txt
exit 0
else 
echo "Input is missing or is in an incorrect format."
exit 1
fi
;;
        "ProteinGBK")
echo "Type in the name of the GBK file"
read -p "GBK file : " gbk
if [ -r $gbk ]; then
echo "Processing..."
python ../../../../opt/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

cat << EOF > pyscript.py
#!/usr/bin/python

import fileinput

for line in fileinput.input():
	tokens = line.split("\t")
	name = tokens[0]
	seq = tokens[1]

	print ">%s\n%s" % (name, seq)

EOF

chmod +x pyscript.py

./pyscript.py Genomebaby.txt > Genome.txt

rm pyscript.py

. ../../../../opt/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. ../../../../opt/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
sed 's,__,%,g' Genomev3.fa | sed 's,_,+,g' | sed 's,+,_1,g' | sed 's,%,_,g' > TiddyGenome.fa
cat TiddyGenome.fa |\
awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' |\
tr "_" "\t" |\
awk -F '   '  '{printf("%s\t%d\n",$0,length($3));}' |\
sort -t '	' -k1,1 -k4,4nr |\
sort -t '	' -k1,1 -u -s |\
sed 's/	/_/' |\
cut -f 1,2 |\
tr "\t" "\n"  |\
fold -w 80 > LongestProtein.fa
grep '^>' LongestProtein.fa | wc -l | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Genes"}' > GeneCount.txt
grep -v '^>' LongestProtein.fa | wc -c | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Aminoacids"}' > AminoacidCount.txt
cat GeneCount.txt AminoacidCount.txt > ProteomeSize.txt
echo "Process completed; please verify the outputs."
echo "Proteome size can be found in 'ProteomeSize.txt'"
rm Genome.txt Genome.fa Genomev2.fa Genomev3.fa TiddyGenome.fa GeneCount.txt AminoacidCount.txt Genomebaby.txt
exit 0 
else 
echo "Input is missing or is in an incorrect format."
exit 1
fi
;;
      *) echo invalid option;; 
esac
done
;;
  "SNPCalling")
            ;;        
        "Help")
            echo "
            ---------------------------Help Manual------------------------------
            Note: Try to run this script in a screen!

	    I. Gene Extraction
              ****************************************************************
              This feature makes use of NCBI's Genbank annotated format files
              (GBKs) to produce a multifasta file containing all ORF's
              belonging to protein coding genes. 
              ****************************************************************
             -1) Collect genbank files belonging to an organism's chromosomes
                 (if assigned); if unassigned, use the single GBK available 
                 (Un). Concatenate all files into a single file.    
             -2) Run MultiTool.sh w/ Extract option. This will output 
                'LongestGenome.fa', containing all ORF's for all protein coding
                 genes.  

            II. Gene Copy Number Analysis
             -1) Run this script w/ Extract option.
             -2) Extract genes of interest from output (LongestGenome.fa).
             -3) Parse 'LongestGenome.fa' to extract desired genes.
             -4) Map reads to individual genes, and then to all ORF's using 
                 GeneCopies option. You should end up with a coverage result 
                 for the individual gene in question, and for all of the ORFs 
                 (Ref:LongestGenome.fa).
             -5) Divide the coverage value belonging to each individual gene 
                 between the coverage value for all ORFs (gene/AllORFs). This
                 is an estimate value for the copy number of that gene.

            III. Find Gene
              This feature extracts a specific gene header and sequence from a
              mutifasta file.
             -1) Run this script w/ Extract option.
             -2) When prompted, use the fasta output from the Extract option as
                 input; specify the gene of interest.

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

            VI. SNP Calling

           VII. PSMC Analysis

            File Naming Protocol (Examples):
               Reference: TP53.fa
               BuildName: TP53bowtie2
               SortedBamName(1rst): TP53.sorted
               SortedBamName(2nd): TP53.sorted.bam
               ResultsName: TP53Cov.txt
               etc.

                              More will be added soon!
            
            Authors: Charles J. Sanfiorenzo Cruz 
                     Jenelys Ruiz Ortiz
                  " 
            ;;
        "Quit")
            break
            ;;
        *) echo invalid option;;
    esac
done

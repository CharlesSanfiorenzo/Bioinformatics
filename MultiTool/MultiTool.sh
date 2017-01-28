#!/bin/bash

#Dependencies: Bowtie2, SAMtools, BBTools, gb2tab, standard unix, python-2.7+, perl 

PS3='Enter your choice (use number): '
options=("Extract" "GeneCopies" "FindGene" "ReadQuality" "ReadFilter" "ProteomeSize" "BaseCounter" "SNPCalling" "MutationLoad" "ToolBox" "Help" "Quit")
alignment=("MultipleAlignments" "SingleAlignment")
suboptions=("SingleRead" "ReadPair")
suboptionspsize=("RnaGBK" "ProteinGBK")
genesoptions=("Exclusive" "Inclusive")
toolbox=("FileDeleter" "Concatenator")
extractoptions=("All Variants" "Longest Variant")
basecountopt=("SingleSeq" "MultiSeq")
alnoptions=("Bowtie2" "Bwa-mem")

#Path to most programs. Modify if needed.
dirpath="../../../../opt"

select opt in "${options[@]}"
do
    case $opt in
        "Extract")
select extropt in "${extractoptions[@]}"
do
    case $extropt in

   "All Variants")
echo "Type in the name of the GBK file"
read -p "GBK file : " gbk
read -p "Output : " output
echo "Processing..."
if [ -r $gbk ]; then
python $dirpath/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

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

. $dirpath/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. $dirpath/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
sed 's,__,%,g' Genomev3.fa | sed 's,_,+,g' | sed 's,+,_1,g' | sed 's,%,_,g' > $output
echo "Process completed; please verify the outputs. Output can be found in 
      '$output'."
rm Genome.txt Genome.fa Genomev2.fa Genomev3.fa Genomebaby.txt
exit 0
else 
echo "Input is missing or is in an incorrect format."
exit 1
fi
;;

 "Longest Variant")
echo "Type in the name of the GBK file"
read -p "GBK file : " gbk
read -p "Output : " output
echo "Processing..."
if [ -r $gbk ]; then
python $dirpath/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

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

. $dirpath/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. $dirpath/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
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
fold -w 80 > temp.out
sed -e  '/^>/s/_.*//' temp.out > $output
echo "Process completed; please verify the outputs. Output can be found in 
      '$output'."
rm Genome.txt Genome.fa Genomev2.fa Genomev3.fa TiddyGenome.fa Genomebaby.txt temp.out
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
        "GeneCopies")
clear
select align in "${alignment[@]}"
do
    case $align in
        "MultipleAlignments")

select subopt in "${suboptions[@]}"
do
    case $subopt in
        "SingleRead")
echo "Number of reference files in directories and/or subdirectories:"
find . -mindepth 1 -type f -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn" -exec printf x \; | wc -c
echo "Full list: "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
echo "Enter number of alignments to be made"
read -p "AlignmentNumber : " number
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fastq"
read -p "Read : " reads

COUNTER=1


	while [ $COUNTER -le $number ]; do
echo -e "\n"
echo "$(tput setaf 6)Alignment $COUNTER"
echo "$(tput setaf 0)Select your reference (.fa) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref_$COUNTER
echo "Select a name for your build and cov results."
read -p "BuildName : " build_$COUNTER
read -p "ResultsName : " CovResults_$COUNTER
echo "Note: If previous values for coverage resulted null, choose 'Yes'"
read -r -p "Do you wish to allow mixed mapping? [y/N] : " response_$COUNTER

let COUNTER=COUNTER+1

	done

COUNTER=1

	while [ $COUNTER -le $number ]; do

cmd="sort"
eval ref=\$ref_$COUNTER
eval build=\$build_$COUNTER
eval CovResults=\$CovResults_$COUNTER
eval response=\$response_$COUNTER

if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $ref ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $reads | bowtie2 -x $build -D 0 -R 0 -N 0 --no-1mm-upfront --no-discordant $reads --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo "Process complete. Please verify '$CovResults'."
exit 0
elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $ref ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $reads | bowtie2 -x $build -D 0 -R 0 -N 0 --no-1mm-upfront --no-mixed --no-discordant $reads --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam $build.sorted.bam
echo "Process complete. Please verify '$CovResults'."
else 
echo "Outputs missing and/or truncated, please verify inputs"
exit 1
fi

let COUNTER=COUNTER+1
done

 ;;
        "ReadPair")

echo "Number of reference files in directories and/or subdirectories:"
find . -mindepth 1 -type f -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn" -exec printf x \; | wc -c
echo "Full list: "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
echo "Enter number of alignments to be made"
read -p "AlignmentNumber : " number
COUNTER=1


	while [ $COUNTER -le $number ]; do
echo -e "\n"
echo "$(tput setaf 6)Alignment $COUNTER"
echo "$(tput setaf 0)Select your reference (.fa) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref_$COUNTER
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fastq"
read -p "Read_1 : " read_1_$COUNTER
read -p "Read_2 : " read_2_$COUNTER
echo "Select a name for your build and cov results."
read -p "BuildName : " build_$COUNTER
read -p "ResultsName : " CovResults_$COUNTER
echo "Note: If previous values for coverage resulted null, choose 'Yes'"
read -r -p "Do you wish to allow mixed mapping? [y/N] : " response_$COUNTER

let COUNTER=COUNTER+1

	done

COUNTER=1

	while [ $COUNTER -le $number ]; do

cmd="sort"
eval ref=\$ref_$COUNTER
eval read_1=\$read_1_$COUNTER
eval read_2=\$read_2_$COUNTER
eval build=\$build_$COUNTER
eval CovResults=\$CovResults_$COUNTER
eval response=\$response_$COUNTER

if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $ref ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $read_1 $read_2 | bowtie2 -x $build -D 0 -R 0 -N 0 --no-1mm-upfront --no-discordant -1 $read_1 -2 $read_2 --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo "Process complete. Please verify '$CovResults'."
exit 0
elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $ref ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $read_1 $read_2 | bowtie2 -x $build -D 0 -R 0 -N 0 --no-1mm-upfront --no-mixed --no-discordant -1 $read_1 -2 $read_2 --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam $build.sorted.bam
echo "Process complete. Please verify '$CovResults'."
else 
echo "Outputs missing and/or truncated, please verify inputs"
exit 1
fi

let COUNTER=COUNTER+1
done

;;

*) echo invalid option;;

esac
done
;;

  "SingleAlignment")

select subopt in "${suboptions[@]}"
do
    case $subopt in
        "SingleRead")
clear
echo "Select your reference (.fa) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fastq"
read -p "Read : " reads
echo "Select a name for your build and cov results."
read -p "BuildName : " build
read -p "ResultsName : " CovResults
echo "Note: If previous values for coverage resulted null, choose 'Yes'"
read -r -p "Do you wish to allow mixed mapping? [y/N] : " response
cmd="sort"
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $ref ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $reads | bowtie2 -x $build -D 0 -R 0 -N 0 --no-1mm-upfront --no-discordant -1 $reads --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo "Process complete. Please verify '$CovResults'."
exit 0
elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $ref ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $reads | bowtie2 -x $build -D 0 -R 0 -N 0 --no-1mm-upfront --no-mixed --no-discordant $reads --threads 16 | samtools view -bS - > $build.bam
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
clear
echo "Select your reference (.fa) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fastq"
read -p "Read_1 : " read_1
read -p "Read_2 : " read_2
echo "Select a name for your build and cov results."
read -p "BuildName : " build
read -p "ResultsName : " CovResults
echo "Note: If previous values for coverage resulted null, choose 'Yes'"
read -r -p "Do you wish to allow mixed mapping? [y/N] : " response
cmd="sort"
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $ref ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $read_1 $read_2 | bowtie2 -x $build -D 0 -R 0 -N 0 --no-1mm-upfront --no-discordant -1 $read_1 -2 $read_2 --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo "Process complete. Please verify '$CovResults'."
exit 0
elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $ref ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $read_1 $read_2 | bowtie2 -x $build -D 0 -R 0 -N 0 --no-1mm-upfront --no-mixed --no-discordant -1 $read_1 -2 $read_2 --threads 16 | samtools view -bS - > $build.bam
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

              *) echo invalid option;;
esac
done
           ;;
    "FindGene")
select genesopt in "${genesoptions[@]}"
do
    case $genesopt in
        "Exclusive")
echo "Select your multifasta and desired gene to be extracted"
read -p "FastaFile : " multifasta
read -p "GeneName : " genename
clear
echo "Processing..."
if [ -r $multifasta ]; then
echo $multifasta $genename | awk 'BEGIN {RS=">"} /'^$genename'/ {print ">"$0}' $multifasta | awk -v RS='>' 'NR>1 {gsub("\n", ";", $0); sub(";$", "", $0); print ">"$0}' | head -n 1 | tr ';' '\n' | sed '/^$/d' > $genename.fa
echo "Process complete. Output can be found in '$genename.fa'."
exit 0
else
echo "Input is missing and/or truncated"
exit 1
fi
           ;;
        "Inclusive")
echo "Select your multifasta and desired gene to be extracted"
read -p "FastaFile : " multifasta
read -p "GeneName : " genename
clear
echo "Processing..."
if [ -r $multifasta ]; then
echo $multifasta $genename | awk 'BEGIN {RS=">"} /'$genename'/ {print ">"$0}' $multifasta | sed '/^$/d' > $genename.fa
echo "Process complete. Output can be found in '$genename.fa'."
exit 0
else
echo "Input is missing and/or truncated"
exit 1
fi
;;
*) echo invalid option;;
        
esac
done
           ;;

 "ReadQuality")
select subopt in "${suboptions[@]}"
do
    case $subopt in
        "SingleRead")
clear
echo "Select your reference (.fa) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fastq"
read -p "Read : " reads
echo "Select a name for your build, coverage & vcf results."
read -p "BuildName : " build
read -p "CovResultsName : " CovResults
read -p "VcfResultsName : " VcfResults
read -r -p "Do you wish to keep files for downstream analysis? [y/N] : " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $ref ] ; then
cmd="sort"
echo $ref $build | bowtie2-build $ref $build
echo $build $reads | bowtie2 -x $build $reads --threads 16 -S $build.sam > AlignmentStats.txt
echo $build | samtools view -bS $build.sam -@ 16 > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > CoverageStats.txt 
cat AlignmentStats.txt CoverageStats.txt > $CovResults
rm AlignmentStats.txt CoverageStats.txt
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.vcf
echo $build $VcfResults | bcftools stats $build.vcf > $VcfResults
exit 0
elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $ref ] ; then
cmd="sort"
echo $ref $build | bowtie2-build $ref $build
echo $build $reads | bowtie2 -x $build $reads --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $build | rm $build.bam
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.vcf
echo $build $VcfResults | bcftools stats $build.vcf > $VcfResults
echo $build | rm $build.sam $build.bam $build.sorted.bam
exit 0
else 
echo "Outputs missing and/or truncated, please verify inputs"
exit 1
fi       
            ;;
        "ReadPair")
clear
echo "Select your reference (.fa) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fastq"
read -p "Read_1 : " read_1
read -p "Read_2 : " read_2
echo "Select a name for your build, coverage & vcf results."
read -p "BuildName : " build
read -p "CovResultsName : " CovResults
read -p "VcfResultsName : " VcfResults
read -r -p "Do you wish to keep files for downstream analysis? [y/N] : " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $ref ] ; then
cmd="sort"
echo $ref $build | bowtie2-build $ref $build
echo $build $read_1 $read_2 | bowtie2 -x $build -1 $read_1 -2 $read_2 --threads 16 -S $build.sam > AlignmentStats.txt
echo $build | samtools view -bS $build.sam -@ 16 > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > CoverageStats.txt
cat AlignmentStats.txt CoverageStats.txt > $CovResults
rm AlignmentStats.txt CoverageStats.txt
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.vcf
echo $build $VcfResults | bcftools stats $build.vcf > $VcfResults
exit 0
elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $ref ] ; then
cmd="sort"
echo $ref $build | bowtie2-build $ref $build
echo $build $read_1 $read_2 | bowtie2 -x $build -1 $read_1 -2 $read_2 --threads 16 | samtools view -bS - > $build.bam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build $CovResults | samtools depth $build.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $CovResults 
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.vcf
echo $build $VcfResults | bcftools stats $build.vcf > $VcfResults
echo $build | rm $build.sam $build.bam $build.sorted.bam
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
 "ReadFilter")
select subopt in "${suboptions[@]}"
do
    case $subopt in
        "SingleRead")
ls -L | find . -name "*.fastq"
read -p "Read : " read
read -p "Minimum Quality Left: " minqual_left
read -p "Minimum Quality Right: " minqual_right
read -p "Minimum Length : " minlength
read -p "Maximum Length : " maxlength
read -p "Trim to : " trim
if [ -r $read ] ; then
perl $dirpath/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $read -trim_qual_right $minqual_right -trim_qual_lef $minqual_left -min_len $minlength -max_len $maxlength -trim_to_len $trim -out_good $read.filtered
exit 0
else 
echo "$read is missing or truncated"
exit 1
fi
;;
	 "ReadPair")
ls -L | find . -name "*.fastq"
read -p "Read_1 : " read_1
read -p "Read_2 : " read_2
read -p "Minimum Quality Left: " minqual_left
read -p "Minimum Quality Right: " minqual_right
read -p "Minimum Length : " minlength
read -p "Maximum Length : " maxlength
read -p "Trim to : " trim
if [ -r $read_1 ] && [ -r $read_2 ] ; then
perl $dirpath/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $read_1 -fastq2 $read_2 -trim_qual_right $minqual_right -trim_qual_lef $minqual_left -min_len $minlength -max_len $maxlength -trim_to_len $trim -out_good filtered
exit 0
else 
echo "$read_1 and $read_2 are missing or truncated"
exit 1
fi
	   ;;
	   *) echo invalid option
	   ;;
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
read -r -p "Do you wish to evaluate longest variant per gene only? [y/N] : " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $gbk ] ; then
echo "Processing..."
python $dirpath/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

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



. $dirpath/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. $dirpath/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
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
cat GeneCount.txt AminoacidCount.txt > LongestProteomeSize.txt
echo "Process completed; please verify the outputs."
echo "Proteome size can be found in 'LongestProteomeSize.txt'"
rm Genome.txt Genome.fa Genomev2.fa Genomev3.fa TiddyGenome.fa GeneCount.txt AminoacidCount.txt Genomebaby.txt LongestProtein.fa
exit 0
elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $gbk ] ; then
echo "Processing..."
python $dirpath/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

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



. $dirpath/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. $dirpath/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
sed 's,__,%,g' Genomev3.fa | sed 's,_,+,g' | sed 's,+,_1,g' | sed 's,%,_,g' > TiddyGenome.fa

grep '^>' TiddyGenome.fa | wc -l | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Genes"}' > GeneCount.txt
grep -v '^>' TiddyGenome.fa | wc -c | awk '{print $1/3}' |  awk '{printf "%2.1f\n",$0}' | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Aminoacids"}' > AminoacidCount.txt
cat GeneCount.txt AminoacidCount.txt > ProteomeSize.txt
echo "Process completed; please verify the outputs."
echo "Proteome size can be found in 'ProteomeSize.txt'"
rm Genome.txt Genome.fa Genomev2.fa Genomev3.fa GeneCount.txt AminoacidCount.txt Genomebaby.txt TiddyGenome.fa
exit 0

else 
echo "Input is missing or is in an incorrect format."
exit 1
fi
;;
        "ProteinGBK")
echo "Type in the name of the GBK file"
read -p "GBK file : " gbk
read -r -p "Do you wish to evaluate longest variant per gene only? [y/N] : " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $gbk ] ; then
echo "Processing..."
python $dirpath/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

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

. $dirpath/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. $dirpath/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
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
cat GeneCount.txt AminoacidCount.txt > LongestProteomeSize.txt
echo "Process completed; please verify the outputs."
echo "Proteome size can be found in 'LongestProteomeSize.txt'"
rm Genome.txt Genome.fa Genomev2.fa Genomev3.fa TiddyGenome.fa GeneCount.txt AminoacidCount.txt Genomebaby.txt LongestProtein.fa
exit 0 
elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $gbk ] ; then
echo "Processing..."
python $dirpath/gb2tab-1.2.1/gb2tab.py -f CDS --genename --block-chars=[E] < $gbk | cut -f1,2 > Genomebaby.txt

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

. $dirpath/bbmap/reformat.sh in=Genome.txt out=Genome.fa ow
perl -p -e 's/^(>.*)$/$1\_/g' Genome.fa > Genomev2.fa
. $dirpath/bbmap/reformat.sh in=Genomev2.fa out=Genomev3.fa uniquenames ow
sed 's,__,%,g' Genomev3.fa | sed 's,_,+,g' | sed 's,+,_1,g' | sed 's,%,_,g' > TiddyGenome.fa
grep '^>' TiddyGenome.fa | wc -l | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Genes"}' > GeneCount.txt
grep -v '^>' TiddyGenome.fa | wc -c | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Aminoacids"}' > AminoacidCount.txt
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
  "BaseCounter")
clear
echo "Type in the name of the fasta file"
read -p "Fasta file : " fasta
if [ -r $fasta ]; then
grep -v '^>' $fasta | wc -c | sed -re ' :restart ; s/([0-9])([0-9]{3})($|[^0-9])/\1,\2\3/ ; t restart ' | awk '{print $1" Bases"}'
echo "Process completed"
exit 0
else 
echo "Input is missing or is in an incorrect format."
exit 1
fi
;;
  "SNPCalling")
select alnopt in "${alnoptions[@]}"
do
    case $alnopt in
        "Bowtie2")
select subopt in "${suboptions[@]}"
do
    case $subopt in

	"SingleRead")
echo "Select your reference (.fna) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref
ls -L | find . -name "*.fastq"
read -p "Read : " read
clear
echo "Select a name for your build, sortedbam, and vcf."
read -p "BuildName : " build
read -r -p "Do you wish to run GATK as well? [y/N] : " response
cmd="sort"
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $ref ] && [ -r $read ] ; then

#Bowtie2 Alignment
echo $ref $build | bowtie2-build $ref $build
echo $build $read | bowtie2 -x $build $read --threads 16 -S $build.sam > $build.aln.txt
echo $build | samtools view -bS $build.sam -@ 16 > $build.bam

#GATK VCF:
echo $ref | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar CreateSequenceDictionary R= $ref O= $ref.dict
echo $build | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar AddOrReplaceReadGroups  INPUT= $build.sam OUTPUT= $build.fixed.sam RGID=1  RGLB= library1  RGPL=illumina RGPU=1  RGSM=sample1 SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
echo $build | rm $build.sam
echo $build | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar SortSam I= $build.fixed.sam  O= $build.fixed.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
echo $build | rm $build.fixed.sam
echo $ref $build | java -Djava.io.tmpdir=../../tmp -jar $dirpath/GenomeAnalysisTK.jar -nct 16 -T HaplotypeCaller -R $ref -I $build.fixed.sorted.bam -stand_call_conf 30 -stand_emit_conf 10 -o $build.gatk.vcf
echo $build | rm $build.fixed.sorted.bam
#Filter
echo $build | vcftools --vcf $build.gatk.vcf --recode --recode-INFO-all --out $build.gatk --minQ 20.00
echo $build | bcftools stats $build.gatk.recode.vcf > $build.filteredbowtiegatkvcf.txt
echo $build | rm $build.gatk.vcf

#Samtools VCF:
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build | rm $build.bam
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.samtools.vcf
echo $build | rm $build.sorted.bam
#Filter
echo $build | vcftools --vcf $build.samtools.vcf --recode --recode-INFO-all --out $build.samtools --minQ 20.00
echo $build | bcftools stats $build.samtools.recode.vcf > $build.filteredbowtiesamtoolsvcf.txt
echo $build | rm $build.samtools.vcf
exit 0

elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $ref ] && [ -r $read ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $read | bowtie2 -x $build $read --threads 16 -S $build.sam > $build.aln.txt
echo $build | samtools view -bS $build.sam > $build.bam
echo $build | rm $build.sam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build | rm $build.bam
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.samtools.vcf
echo $build | rm $build.sorted.bam 
echo $build | vcftools --vcf $build.samtools.vcf --recode --recode-INFO-all --out $build.samtools --minQ 20.00
echo $build | bcftools stats $build.samtools.recode.vcf > $build.filteredbowtiesamtoolsvcf.txt
echo $build | rm $build.samtools.vcf
exit 0

else 
echo "Outputs missing and/or truncated, please verify inputs"
exit 1
fi

;;

	"ReadPair")
echo "Select your reference (.fna) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref
ls -L | find . -name "*.fastq"
read -p "Read_1 : " read_1
read -p "Read_2 : " read_2
clear
echo "Select a name for your build, sortedbam, and vcf."
read -p "BuildName : " build
read -r -p "Do you wish to run GATK as well? [y/N] : " response
cmd="sort"
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] && [ -r $ref ] && [ -r $read_1 ] && [ -r $read_2 ] ; then

#Bowtie2 Alignment
echo $ref $build | bowtie2-build $ref $build
echo $build $read_1 $read_2 | bowtie2 -x $build -1 $read_1 -2 $read_2 --threads 16 -S $build.sam > $build.aln.txt
echo $build | samtools view -bS $build.sam -@ 16 > $build.bam

#GATK VCF:
echo $ref | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar CreateSequenceDictionary R= $ref O= $ref.dict
echo $build | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar AddOrReplaceReadGroups  INPUT= $build.sam OUTPUT= $build.fixed.sam RGID=1  RGLB= library1  RGPL=illumina RGPU=1  RGSM=sample1 SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
echo $build | rm $build.sam
echo $build | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar SortSam I= $build.fixed.sam  O= $build.fixed.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
echo $build | rm $build.fixed.sam
echo $ref $build | java -Djava.io.tmpdir=../../tmp -jar $dirpath/GenomeAnalysisTK.jar -nct 16 -T HaplotypeCaller -R $ref -I $build.fixed.sorted.bam -stand_call_conf 30 -stand_emit_conf 10 -o $build.gatk.vcf
echo $build | rm $build.fixed.sorted.bam
#Filter
echo $build | vcftools --vcf $build.gatk.vcf --recode --recode-INFO-all --out $build.gatk --minQ 20.00
echo $build | bcftools stats $build.gatk.recode.vcf > $build.filteredbowtiegatkvcf.txt
echo $build | rm $build.gatk.vcf

#Samtools VCF:
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build | rm $build.bam
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.samtools.vcf
echo $build | rm $build.sorted.bam
#Filter
echo $build | vcftools --vcf $build.samtools.vcf --recode --recode-INFO-all --out $build.samtools --minQ 20.00
echo $build | bcftools stats $build.samtools.recode.vcf > $build.filteredbowtiesamtoolsvcf.txt
echo $build | rm $build.samtools.vcf
exit 0

elif [[ $response =~ ^([nN][oO]|[nN])$ ]] && [ -r $ref ] && [ -r $read ] ; then
echo $ref $build | bowtie2-build $ref $build
echo $build $read_1 $read_2 | bowtie2 -x $build -1 $read_1 -2 $read_2 --threads 16 -S $build.sam > $build.aln.txt
echo $build | samtools view -bS $build.sam -@ 16 > $build.bam
echo $build | rm $build.sam
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build | rm $build.bam
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.samtools.vcf
echo $build | rm $build.sorted.bam 
echo $build | vcftools --vcf $build.samtools.vcf --recode --recode-INFO-all --out $build.samtools --minQ 20.00
echo $build | bcftools stats $build.samtools.recode.vcf > $build.filteredbowtiesamtoolsvcf.txt
echo $build | rm $build.samtools.vcf
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
	"Bwa-mem")

echo "Select your reference (.fna) and read (.fastq) files; include path."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls -L | find . -name "*.fa" -o -name "*.faa" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn"
read -p "Reference : " ref
ls -L | find . -name "*.fastq"
echo "For a single read alignment, leave 'Read_2' empty and hit enter"
read -p "Read_1 : " read_1
read -p "Read_2 : " read_2
clear
echo "Select a name for your build, sortedbam, and vcf."
read -p "BuildName : " build
cmd="sort"
if [ -r $ref ] ; then

#Bwa mem alignment
echo $ref | bwa index -a $ref
echo $ref $read_1 $read_2 $build | bwa mem -M $ref $read_1 $read_2 > $build.sam
echo $build | samtools view -bS $build.sam > $build.bam

#GATK VCF:
echo $ref | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar CreateSequenceDictionary R= $ref O= $ref.dict
echo $build | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar AddOrReplaceReadGroups  INPUT= $build.sam OUTPUT= $build.fixed.sam RGID=1  RGLB= library1  RGPL=illumina RGPU=1  RGSM=sample1 SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
echo $build | rm $build.sam
echo $build | java -Djava.io.tmpdir=../../tmp -XX:ParallelGCThreads=16 -jar $dirpath/picard-tools-2.0.1/picard.jar SortSam I= $build.fixed.sam  O= $build.fixed.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
echo $build | rm $build.fixed.sam
echo $ref $build | java -Djava.io.tmpdir=../../tmp -jar $dirpath/GenomeAnalysisTK.jar -nct 16 -T HaplotypeCaller -R $ref -I $build.fixed.sorted.bam -stand_call_conf 30 -stand_emit_conf 10 -o $build.gatk.vcf
echo $build | rm $build.fixed.sorted.bam

#Filtering based on Quality
echo $build | vcftools --vcf $build.gatk.vcf --recode --recode-INFO-all --out $build.bwagatk --minQ 20.00
echo $build | bcftools stats $build.bwagatk.recode.vcf > $build.filteredbwagatkvcf.txt
echo $build | rm $build.gatk.vcf

#Samtools VCF:
echo $cmd $build | samtools $cmd $build.bam $build.sorted -@ 16
echo $build | rm $build.bam
echo $ref $build | samtools mpileup -uf $ref $build.sorted.bam | bcftools call -c > $build.bwasamtools.vcf
echo $build | rm $build.sorted.bam

#Filter
echo $build | vcftools --vcf $build.bwasamtools.vcf --recode --recode-INFO-all --out $build.bwasamtools --minQ 20.00
echo $build | bcftools stats $build.bwasamtools.recode.vcf > $build.filteredbwasamtoolsvcf.txt
echo $build | rm $build.bwasamtools.vcf

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

  "MutationLoad")

echo "Select your canonical gene references (.fa), SnpEff library, and vcf."
echo "Note: You may choose from the following files in this directory and/or 
subdirectories "
ls ../../snpEff/data/ 
read -p "SnpEffLibrary : " lib
ls -L | find . -name "*.vcf"
read -p "Vcf : " vcf

if [ -r $vcf ] ; then
#Creates the database for the genome to be used
echo $lib | java -Xmx20g -Djava.io.tmpdir=../../tmp -jar ../../snpEff/snpEff.jar build -gff3 $lib -c ../../snpEff/snpEff.config -t -v
#Runs SNPEff
echo $lib $vcf | java -Xmx20g -Djava.io.tmpdir=../../tmp -jar ../../snpEff/snpEff.jar -v -csvStats stats.csv $lib $vcf > $vcf.ann.vcf

#Keeps only the 'Effect' column and count number of synonymous and non-synonymous SNPs.
cat $vcf.ann.vcf | cut -f 8 | tr ";" "\n" | grep ^ANN= | cut -f 2 -d = | tr "," "\n" > temp.out
#Synonymous SNPs
grep -e synonymous_variant -e -start_retained -e stop_retained_variant temp.out | wc -l > temp2.out
#Non-Synonymous SNPs
grep -e missense_variant -e stop_gained -e stop_lost -e initiator_codon_variant -e start_lost temp.out | wc -l > temp3.out
#High impact
grep -e HIGH temp.out | wc -l > temp4.out
#Low impact
grep -e LOW temp.out | wc -l > temp5.out
#Moderate impact
grep -e MODERATE temp.out | wc -l > temp6.out
#Protein Coding Count
grep -e protein_coding temp.out | wc -l > temp7.out


#Then adapt for upstream analysis (python config parser module
sed -i '1s/^/non-synonymous_snp : /' temp3.out
sed -i '1s/^/synonymous_snp : /' temp2.out
sed -i '1s/^/HighImpact : /' temp4.out
sed -i '1s/^/LowImpact : /' temp5.out
sed -i '1s/^/ModerateImpact : /' temp6.out
sed -i '1s/^/Protein coding SNPs : /' temp7.out
rm temp.out 
#Counts SNP syn and non-syn sites, parses SNPEff output, and concatenates them (for pyhton config parser)
perl ../../KsKa.pl > temp.out

cat temp.out temp3.out temp2.out temp4.out temp5.out temp6.out temp7.out > SNPStats.txt
rm temp.out temp2.out temp3.out temp4.out temp5.out temp6.out temp7.out
#Calculate pKa/Ks using Nei-Gojobori
cat << EOF > pyscript.py
#!/usr/bin/python
#This short script uses the output values of KaKs.pl & SnpEff to calculate mutational load using Nei-Gojobori: pKa/Ks = [-3/4ln(1-4pn/3)] / [-3/4ln(1-4ps/3)], where ps = syn SNPs / syn sites and pn = nonsyn SNPs / nonsyn sites
from math import log #If for some reason you need to calculate the logarithm of a negative number, import cmath instead.
import ConfigParser

config = ConfigParser.RawConfigParser()
config.read("SNPStats.txt")
nonSyn_site = float(config.get("myvars", "non-synonymous_number"))
Syn_site = float(config.get("myvars", "synonymous_number"))
nonSyn_SNP = float(config.get("myvars", "non-synonymous_snp"))
Syn_SNP = float(config.get("myvars", "synonymous_snp"))
HighImp = int(config.get("myvars", "HighImpact"))
LowImp = int(config.get("myvars", "LowImpact"))
ModImp = int(config.get("myvars", "ModerateImpact"))

pn = nonSyn_SNP/nonSyn_site
ps = Syn_SNP/Syn_site

print "The pKs/Ks ratio for this organism is:", (-3/4*log(1-(4*pn)/3))/(-3/4*log(1-(4*ps)/3))
print "High Impact Percentage:", (HighImp*100/(HighImp+LowImp+ModImp))
print "Moderate Impact Percentage:", (ModImp*100/(HighImp+LowImp+ModImp))
print "Low Impact Percentage:", (LowImp*100/(HighImp+LowImp+ModImp))
print "Impact average:", (HighImp*3/((HighImp*3)+(LowImp)+(ModImp*2)))

EOF

chmod +x pyscript.py
python pyscript.py > pKaKs.txt
rm pyscript.py
echo "Process complete, please verify 'pKaKs.txt'"
exit 0
else 
echo "Input is missing and/or truncated"
exit 1
fi

            ;;

  "ToolBox")
clear
select tools in "${toolbox[@]}"
do
    case $tools in
  "FileDeleter")
echo "This feature will delete all files within a directory with a given 
extension"
echo "Specify extension; include dot. Example: '.fa'"
read -p "Extension : " ext
read -r -p "This will delete all files with '$ext'. Are you sure? [y/N] : " response2

if [[ $response2 =~ ^([yY][eE][sS]|[yY])$ ]] ; then
echo "..."
echo $ext | find . -name "*$ext" -exec rm -rf {} \;
echo "All files with '$ext' extension have been deleted. 
      Verify directory content."
exit 0
elif [[ $response2 =~ ^([nN][oO]|[nN])$ ]] ; then
echo "Exiting..."
exit 1

else
echo "invalid option"

fi
;;
  "Concatenator")
echo "This feature will concatenate all files within a directory with a given 
extension"
echo "Specify extension; include dot. Example: '.fa'"
read -p "Extension : " ext
read -p "Output : " catname

echo $ext | find . -name "*$ext" -exec cat {} \; > $catname

echo "All files with '$ext' extension have been concatenated. 
      Verify '$catname'."
exit 0
;;
*) echo invalid option;;


esac
done
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
            
           -Authors: 
            Charles J. Sanfiorenzo Cruz 
            Jenelys Ruiz Ortiz

           -Collaborators:

                                                                 Version 1.4 
                  " 
            ;;
        "Quit")
            break
            ;;
        *) echo invalid option;;
    esac
done

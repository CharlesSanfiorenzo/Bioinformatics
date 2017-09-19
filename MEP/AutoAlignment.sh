##############Automated pipeline for integration of datasets.###################
######If you have any questions, please contact the author of this script#######
######at csanfior@mit.edu ######################################################

#Find and verify paths for Bowtie2, TopHat2, and Samtools

##Bowtie2
echo "Verifying that Bowtie2 is installed ..."
#exec 2> /dev/null
bowtie2=$(locate -l 1 "bowtie2")
if [ ! -d $bowtie2 ] ; then                                                   #1
#exec 2> /dev/tty
echo "Could not find a valid Bowtie2 install path. Please specify one."
read -p "Path : " bowtie2
if [ ! -d $bowtie2 ] ; then                                                   #2
echo "Invalid path. Please download and install the newest version of Bowtie2. Exiting now."
exit 1
fi                                                                            #2

else
echo "Found it!"

##TopHat2
echo "Verifying that TopHat2 is installed ..."
#exec 2> /dev/null
tophat=$(locate -l 1 "tophat-")
if [ ! -d $tophat ] ; then                                                    #3
#exec 2> /dev/tty
echo "Could not find a valid TopHat2 install path. Please specify one."
read -p "Path : " tophat
if [ ! -d $tophat ] ; then                                                    #4
echo "Invalid path. Please download and install the newest version of TopHat2. Exiting now."
exit 1
fi                                                                            #4
fi                                                                            #3

echo "Found it!"

##Samtools
echo "Verifying that Samtools is installed ..."
#exec 2> /dev/null
samtools=$(locate -l 1 "samtools")
if [ ! -d $samtools ] ; then                                                  #5
#exec 2> /dev/tty
echo "Could not find a valid Samtools install path. Please specify one."
read -p "Path : " cufflinks
if [ ! -d $samtools ] ; then                                                  #6
echo "Invalid path. Please download and install the newest version of Samtools. Exiting now."
exit 1
fi                                                                            #6
fi                                                                            #5

echo "Found it!"
echo "Starting protocol..."

#TopHat2 Alignment
ls -L | find . -name "*.gtf"
read -p "Gene annotations: " genes
ls -L | find . -name "*.fastq"
read -p "Reads (e.g. read_1.fastq read_2.fastq): " reads
read -p "Output directory: " output

echo "
      *********************************************************************
      Do you wish to build a genome index?
      Note that if this was done in previous runs, this step can be skipped.
      ********************************************************************* "
loop=true
while $loop = "true" ; do                                               #loopStart
read -p "Response [y/N] : " response 
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] ; then                            #7
ls -L | find . -name "*.fa"
read -p "Genome fasta: " genome
$bowtie2/bowtie2-build $genome Genome #Bowtie2 index
if [ ! -d "Transcriptome" ] ; then                                            #8
$tophat/tophat2 --GTF $genes --transcriptome-index Transcriptome/genes -o $output -p 16 --no-novel-juncs -g 1 Genome $reads #TopHat2 Alignment
else
$tophat/tophat2 --transcriptome-index Transcriptome/genes -o $output -p 16 --no-novel-juncs -g 1 Genome $reads #TopHat2 Alignment
fi                                                                            #8
$samtools/samtools bam2fq $output/unmapped.bam > $output/unmapped.fastq #Convert unaligned bam to fastq
$tophat/tophat2 --transcriptome-index Transcriptome/genes -o $output/Unmapped -p 16 --no-novel-juncs -N 5 --read-edit-dist 5 --read-gap-length 5 --b2-N 1 -g 1 Genome $output/unmapped.fastq #Second TopHat2 Alignment
$samtools/samtools merge $output/Merged.bam $output/accepted_hits.bam $output/Unmapped/accepted_hits.bam #Merge bams of accepted hits
loop="false"

elif [[ $response =~ ^([nN][oO]|[nN])$ ]] ; then
if [ ! -d "Transcriptome" ] ; then                                            #9
$tophat/tophat2 --GTF $genes --transcriptome-index Transcriptome/genes -o $output -p 16 --no-novel-juncs -g 1 Genome $reads #TopHat2 Alignment
else
$tophat/tophat2 --transcriptome-index Transcriptome/genes -o $output -p 16 --no-novel-juncs Genome $reads #TopHat2 Alignment
fi                                                                            #9
$samtools/samtools bam2fq $output/unmapped.bam > $output/unmapped.fastq #Convert unaligned bam to fastq
$tophat/tophat2 --transcriptome-index Transcriptome/genes -o $output/Unmapped -p 16 --no-novel-juncs -N 5 --read-edit-dist 5 --read-gap-length 5 --b2-N 1 -g 1 Genome $output/unmapped.fastq #Second TopHat2 Alignment
$samtools/samtools merge $output/Merged.bam $output/accepted_hits.bam $output/Unmapped/accepted_hits.bam #Merge bams of accepted hits
loop="false"

else
echo "Please provide a valid answer..."

fi                                                                            #7
done                                                                    #loopEnd

fi                                                                            #1

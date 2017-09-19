##############Automated pipeline for integration of datasets.###################
######If you have any questions, please contact the author of this script#######
######at csanfior@mit.edu ######################################################

#Find and verify paths for Bowtie2, TopHat2, and Samtools

##Bowtie2
echo "Verifying that Cufflinks is installed ..."
#exec 2> /dev/null
cufflinks=$(locate -l 1 "cufflinks")
if [ ! -d $cufflinks ] ; then                                                 #1
#exec 2> /dev/tty
echo "Could not find a valid Cufflinks install path. Please specify one."
read -p "Path : " cufflinks
if [ ! -d $cufflinks ] ; then                                                 #2
echo "Invalid path. Please download and install the newest version of Cufflinks. Exiting now."
exit 1
fi                                                                            #2

else

echo "Found it!"
echo "Starting protocol..."

#Cufflinks expression analysis
ls -L | find . -name "*.gtf"
read -p "Gene annotations: " genes
ls -L | find . -name "*.fa"
read -p "Genome Assembly: " assembly
read -p "Control Ribo-Seq TopHat Results Directory : " Ribo1
read -p "Experimental Ribo-Seq TopHat Results Directory : " Ribo2
read -p "Control RNA-Seq TopHat Results Directory : " RNA1
read -p "Experimental RNA-Seq TopHat Results Directory : " RNA2
read -p "Output directory : " output
echo "
      *********************************************************************
      Do you wish to select genes from a list of gene symbols?
      Note that if no gene symbol list is used, all genes will be processed.
      ********************************************************************* "
read -p "Response [y/N] : " response

while [[ $response != ^([yY][eE][sS]|[yY]|[nN][oO]|[nN])$ ]] ; do
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] ; then
read -p "Gene Symbol list : " symbols 
fi
done

$cufflinks/cuffdiff -p 16 -b $assembly -o $output $genes $Ribo1/Merged.bam $Ribo2/Merged.bam $RNA1/Merged.bam $RNA2/Merged.bam #Cufflinks run

awk '{ print $1, $10, $14, $18, $22}' $output/genes.fpkm_tracking | awk '{if (NR!=1) {print}}' > $output/AllGenesFPKM.txt #Getting columns we want ...

if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]] ; then
sed 's,-,_,g' $output/AllGenesFPKM.txt | grep -Fwf $symbols | sed 's,_,-,g' > $output/SelectGenesFPKM.txt #Extract specific genes ^ we need to remove first row
fi 

fi                                                                            #1

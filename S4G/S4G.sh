#! /bin/bash
#Author: Charles Sanfiorenzo
#Files you will need: genome.fna, genome.gff, & genome.gbff. Easily found in NCBI ftp
#Tools you will need: standard GNU tools, Bedtools, Python 2.7 (w/ Pandas)

#Note: Extreme downloads: curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Helicobacter_pylori/assembly_summary.txt' | awk '{FS="\t"} !/^#/ {print $20}' | sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Helicobacter_pylori/latest_assembly_versions/)(GCF/)([0-9]{3}/)([0-9]{3}/)([0-9]{3}/)(GCF_.+)|\1\2\3\4\5\6/\6_genomic.fna.gz|' | sed 's/$/_genomic.fna.gz/' > genomic_file

if [ "$1" == '' ] || [ "$2" == '' ] || [ "$3" == '' ] ; then
	echo 'Usage: ./S4G.sh genome.fna genome.gbff genome.gff
              Follow that specific order!'
	exit 1
elif [ ! -s "$1" ] || [ ! -s "$2" ] || [ ! -s "$3" ] ; then
	echo 'Usage: ./S4G.sh genome.fna genome.gbff genome.gff
              Follow that specific order!'
	echo 'One or more of the files specified are not present in this directory...'
	exit 1
fi
#Emergency fill
gff="$3"
gbff="$2"
fna="$1"

###Intron Length

#Removes genes with ncRNA classification (that do not fit within mRNA, tRNA, snRNA, microRNA, and/or rRNA classifications). Only leaves exon annotations for rest of genes. Removes low quality sequences
grep -v 'ncRNA' $gff | awk -F'\t' '$3!="gene"' | awk -F'\t' '$3!="mRNA"' | awk -F'\t' '$3!="CDS"' | awk -F'\t' '$3!="region"' | awk -F'\t' '$3!="gap"' | awk -F'\t' '$3!="sequence_feature"' | awk -F'\t' '$3!="tRNA"' | awk -F'\t' '$3!="rRNA"' | grep -v 'low-quality sequence' | grep -v 'pseudogene' > tmp.gff
#Converts gffs into bed:
awk -F '\t' '{ print $9,'\t',$4,'\t',$5,'\t',$7 }' tmp.gff | awk -F 'transcript_id=' '{print $NF}' | grep -v '^#' | sort | sed 's,  ,\t,g'  > $gff.noNcRNA_Exons.bed

#Converts gffs into bed for intergenic:
awk -F '\t' '{ print $1,'\t',$4,'\t',$5,'\t',$7 }' tmp.gff | awk -F 'transcript_id=' '{print $NF}' | grep -v '^#' | sort | sed 's,  ,\t,g'  > $gff.noNcRNA_ExonsV2.bed


if [ -s $gff.noNcRNA_Exons.bed ]; then
#Creates python script. Checks for introns, and count overall intron length within genome
cat << EOF > pyscript.py
#!/usr/bin/python
#Author: Charles Sanfiorenzo
#Gets intron length from a bed file
import pandas as pd
import sys, getopt
import os.path

def main(argv) :
   global inputFile
   try:
      opts, args = getopt.getopt(argv,"f:",[])
   except getopt.GetoptError:
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-f"):
         inputFile = arg
   
if __name__ == "__main__":
	main(sys.argv[1:])

df = pd.read_csv(inputFile, sep="\t")
columns = ['accession','start','end','strand']
df.columns = columns
df.dropna(inplace=True)
df.index = range(len(df.index))
#print df #Diagnostic
counter = 1
oldCounter = 0 #oldCounter is previous accession index
intronSum = 0
intronCount = 0
genes = []
intronSums = []
for i in range(1,len(df.index)) :
	if df.ix[counter]['accession'] == df.ix[oldCounter]['accession'] :
		intronSum += int(df.ix[counter]['start']-df.ix[oldCounter]['end'])
		intronSums += [intronSum]
		intronCount += 1 
		#print df.ix[counter]['accession'] #Diagnostic
		if df.ix[counter]['accession'] not in genes :
			genes += [str(df.ix[counter]['accession'])]
	counter += 1
	oldCounter += 1
print 'Total genes with introns:', len(genes)
print 'Total introns:', intronCount
print 'Total intron length:', intronSum
print 'Average intron length:', float(intronSum)/float(intronCount)
EOF

chmod +x pyscript.py
./pyscript.py -f $gff.noNcRNA_Exons.bed
rm pyscript.py
fi

rm tmp.gff
rm $gff.noNcRNA_Exons.bed

###Intergenic Length (feature extract)
awk -F ' ' '/^>/ { print $1; next } 1' $fna > genome.fa
cat genome.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t0\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > genome.bed
bedtools subtract -a genome.bed -b $gff.noNcRNA_ExonsV2.bed > Intergenic.bed
rm $gff.noNcRNA_ExonsV2.bed
bedtools getfasta -fi genome.fa -bed Intergenic.bed -fo Intergenic.fa
rm Intergenic.bed
rm genome.bed

#Gets nucleotide count from Intergenic.fa
inter=$(grep -v '^>' Intergenic.fa | wc -c)
echo "Intergenic length: $inter"
count=$(grep -v '^>' $fna | wc -c)
size=$count
echo "Genome size: $count"
ratio=$(bc <<< "scale=2; $inter/$count")
echo "Intergenic/Genome Ratio: $ratio"
# Removes header lines, linearizes sequences (removes newlines), and counts CG pairs.
count=$(grep -v '^>' $fna | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n//g' | grep -o 'CG' | wc -l)
echo "CpG sites: $count"
count=$(grep -v '^>' $fna | grep -E 'C|G' | wc -c)
ratio=$(bc <<< "scale=3; $count/$size") 
echo "GC Content: $ratio"

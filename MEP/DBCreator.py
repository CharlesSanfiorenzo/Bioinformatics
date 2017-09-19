######Produces final databse for SVM.###########################################
######If you have any questions, please contact the author of this script#######
######at csanfior@mit.edu ######################################################

import sys, getopt
import pandas as pd
import os.path
import re #Used to retain only numbers from Needle alignment scores
from collections import Counter


#Note Finish adding biasT later
########Argument parser
def main(argv):
   global FPKMfile
   global needle
   global fastaFile
   global added
   global add
   global biasT
   global best
   global noVariants
   FPKMfile = ''
   needle = ''
   fastaFile = ''
   required = []
   added = ''
   add = False
   bias = False
   best = False
   noVariants = False
   try:
      opts, args = getopt.getopt(argv,"hf:M:N:a:",["sample-t:","best-var","no-var"])
   except getopt.GetoptError:
      print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRefGene.fpkm.csv> -N <NeedleResults>'
      print 'Note: Type \'-h\' for additional help'
      sys.exit(2)
   for opt, arg in opts:
      if len(sys.argv) == 1 :
	print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRefGene.fpkm.csv> -N <NeedleResults>'
      	print 'Note: Type \'-h\' for additional help'
      	sys.exit(2)
      elif opt == '-h':
         print '''Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRefGene.fpkm.csv> -N <NeedleResults> [optional: -a <additional feature> --sample-t <bias threshold> --best-var]
                  Author: Charles Sanfiorenzo (2017)
		  --------------------------------------------------------------------
                  Running this script will create a database for MEP2's SVM module. The
		  resulting file (MEP2.db) will be a csv containing the following info:
                  GeneSymbol, MotifSequence, AlignmentScore, and TranslationalEfficiency.

		  Optional Arguments: 
		  -a | Uses a csv with collected features (do not specify column name, it will be asked via prompt)

		  For any additional questions, email the author of this script at
		  csanfior@mit.edu'''
         sys.exit()
      elif opt in ("-f"):
         fastaFile = arg
	 required += [arg]
      elif opt in ("-M"):
         FPKMfile = arg
	 required += [arg]
      elif opt in ("-N"):
         needle = arg
	 required += [arg]
      elif opt in ("-a"):
         added = arg
	 add = True
      #elif opt in ("--sample-t"):
        # biasT = arg
	# bias = True

      elif opt in ("--no-var"):
	 noVariants = True

      elif opt in ("--best-var"):
         best = True
      else :
	 print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRefGene.fpkm.csv> -N <NeedleResults>'
         print 'Note: Type \'-h\' for additional help'
	 sys.exit(2)
   args = sys.argv[1:]
   if not args :
	print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRefGene.fpkm.csv> -N <NeedleResults>'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if len(required) == 1 :
	print 'Error: Missing argument.' 
	print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRefGene.fpkm.csv> -N <NeedleResults>'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(fastaFile) == False :
	print 'Error:',fastaFile,'not found' 
	print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRefGene.fpkm.csv> -N <NeedleResults>'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(FPKMfile) == False :
	print 'Error:',FPKMfile,'not found' 
	print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRefGene.fpkm.csv> -N <NeedleResults>'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(needle) == False :
	print 'Error:',needle,'not found' 
	print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRef.fpkm.csv> -N <NeedleResults>'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if add == True and os.path.isfile(added) == False :
	print 'Error:',added,'not found' 
	print 'Usage: DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRef.fpkm.csv> -N <NeedleResults>'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   #if bias == True and biasT.is_numeric() == False :
	#print 'Error:',biasT,'is not a number' 
	#print 'Usage:  DBCreator.py -f <fasta(no ref)> -M <MEP2_NoRef.fpkm.csv> -N <NeedleResults>'
       # print 'Note: Type \'-h\' for additional help'
	#sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])


########

# A nice little file reader
def contentExtractor(files) :
	with open(files) as f:
    		content = f.readlines()
	content = [x.strip('\n') for x in content] 
	return content

content = contentExtractor(fastaFile)
############ Generator 1 (Header and sequence extractor)
sequences = []
headers = []
for element in content :

	if element[0] != '>' :
		sequences += [element]
	else :
		headers += [element.strip('>')]
############


# Load FPKM annotations
df2 = pd.read_csv(FPKMfile, header=0, sep=",", index_col=0, parse_dates=True)
columns=['Gene Symbol', 'Ribo-Seq: Control (Optimal conditions) FPKM', 'Ribo-Seq: Amino acid starvation FPKM', 'RNA-Seq: Control (Optimal conditions) FPKM', 'RNA-Seq: Amino acid starvation FPKM', 'Translational Efficiency']
df2.columns = columns
df2.reset_index(drop=True,inplace=True)

# Create dataframe for database
idx = range(len(sequences))
df = pd.DataFrame(index=idx)
df["Gene Symbol"] = headers
df["Sequence"] = sequences
df.reset_index(drop=True,inplace=True)
content = contentExtractor(needle)

############# Generator 2 (Alignment score extractor) # Recently added
scores = [s for s in content if "Score:" in s]
numOnly = []
for element in scores :
	numOnly += [re.sub("[^0-9|.]", "", element)]
#############

df["Alignment Score"] = numOnly

############# Generator 3 (Translation efficiency extractor)
#Generate number of variants per gene
counterLst = []
sameCounter = 0
counter = 0
while sameCounter <= len(df.index) :
	newID = sameCounter + 1
	if newID < len(df.index) and df.ix[sameCounter]['Gene Symbol'] != df.ix[newID]['Gene Symbol'] :
		counter += 1
		counterLst += [counter]
		counter = 0
		sameCounter += 1
	elif newID < len(df.index) and df.ix[sameCounter]['Gene Symbol'] == df.ix[newID]['Gene Symbol'] :
		counter += 1
		sameCounter += 1
	else :
		sameCounter += 1
		
#Checksum for last element
if headers[-1] == headers[-2] :
	counterLst[-1] += 1
else :
	counterLst += [1]

#Use list of numbers of variants per gene to write Translation Efficiency the correct amount of times
translationEff = []
ribo1FPKM = []
ribo2FPKM = []
RNA1FPKM = []
RNA2FPKM = []
sameCounter = 0
for number in counterLst :
	for i in range(number) :

		if type(df2.ix[:,"Translational Efficiency"][sameCounter]) == str :
			translationEff += [str(df2.ix[:,"Translational Efficiency"][sameCounter])]
			ribo1FPKM += [float(df2.ix[:,"Ribo-Seq: Control (Optimal conditions) FPKM"][sameCounter])]
			ribo2FPKM += [float(df2.ix[:,"Ribo-Seq: Amino acid starvation FPKM"][sameCounter])]
			RNA1FPKM += [float(df2.ix[:,"RNA-Seq: Control (Optimal conditions) FPKM"][sameCounter])]
			RNA2FPKM += [float(df2.ix[:,"RNA-Seq: Amino acid starvation FPKM"][sameCounter])]
			i += 1
		else :
			translationEff += [float(df2.ix[:,"Translational Efficiency"][sameCounter])]
			ribo1FPKM += [float(df2.ix[:,"Ribo-Seq: Control (Optimal conditions) FPKM"][sameCounter])]
			ribo2FPKM += [float(df2.ix[:,"Ribo-Seq: Amino acid starvation FPKM"][sameCounter])]
			RNA1FPKM += [float(df2.ix[:,"RNA-Seq: Control (Optimal conditions) FPKM"][sameCounter])]
			RNA2FPKM += [float(df2.ix[:,"RNA-Seq: Amino acid starvation FPKM"][sameCounter])]
			i += 1
	sameCounter += 1

#############
df["Ribo-Seq: Control (Optimal conditions) FPKM"] = ribo1FPKM
df["Ribo-Seq: Amino acid starvation FPKM"] = ribo2FPKM
df["RNA-Seq: Control (Optimal conditions) FPKM"] = RNA1FPKM
df["RNA-Seq: Amino acid starvation FPKM"] = RNA2FPKM
df["Translational Efficiency"] = translationEff
df.reset_index(drop=True,inplace=True) 



#Only keep genes without variants
if noVariants == True :
	best == False
	df.drop_duplicates(subset='Gene Symbol', inplace=True,keep=False)
	df.reset_index(drop=True,inplace=True)


#Keep only best score variant
#Use list of numbers of variants per gene to write Translation Efficiency the correct amount of times

if best == True :

	masker = df.groupby(['Gene Symbol'])['Alignment Score'].transform(max) == df['Alignment Score']
	

# Feature to be added
if added != '' and add == True :
	content = contentExtractor(added)
	names = str(raw_input("Name of added feature: ")) #Use normal input for Python3
	df[str(names)] = content
	with open("MEP2.settings", "a") as myfile:
		myfile.write("Added Feature = Yes\nFeature = "+names)
	#Makeshift solution for python 2.7
	#fin = open('MEP2.settings', "r" )
	#data_list = fin.readlines()
	#fin.close()
	#del data_list[-1:-9]
	#fout = open('MEP2.settings', "w")
	#fout.writelines(data_list)
	#fout.close()
elif added == '' :
	with open("MEP2.settings", "a") as f:
		f.write("Added Feature : No")

print 'Process complete: Verify MEP2.db'

if best == True :
	df[masker].to_csv('MEP2.db')
else :
	df.to_csv('MEP2.db')

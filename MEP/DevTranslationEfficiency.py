######Calculates translation efficiency and updates database & fasta.###########
######If you have any questions, please contact the author of this script####### #Add if upregulation or downregulation
######at csanfior@mit.edu ######################################################


import sys, getopt
import pandas as pd
import os.path

########Argument parser
def main(argv):
   global FPKMfile
   global fastaFile
   global cutaway
   global regType
   global bias
   regType = 'down'
   FPKMfile = ''
   fastaFile = ''
   cutaway = 33.0
   required = []
   bias = False
   try:
      opts, args = getopt.getopt(argv,"hT:f:r:b")
   except getopt.GetoptError:
      print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
      print 'Note: Type \'-h\' for additional help'
      sys.exit(2)
   for opt, arg in opts:
      if len(sys.argv) == 1 :
	print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
      	print 'Note: Type \'-h\' for additional help'
      	sys.exit(2)
      elif opt == '-h':
         print '''Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa
		  Example: TranslationalEfficiency.py -T SelectGenesFPKM.txt -f SelectGenesMotif.fa -r 33 -R down > MotifsNoRef.fa
                  Author: Charles Sanfiorenzo (2017)
		  --------------------------------------------------------------------
                  Running this script will calculate translation efficiency and produce 
                  a csv with the results (MEP.fpkm.csv) in the working directory. The
		  gene with the most downregulation is used as reference for alignment
		  score determination, and so an updated csv (MEP_NoRefGenes.fpkm.csv) and
		  two fastas (ReferenceMotif.fa & MotifsNoRef.fa) are created for later
		  steps in the protocol. 
                  
                  Optional arguments:
                  -r | establishes the allowed difference in percentage of RNA transcription
                       per gene. Default: 33
                  -R | Regulation type under study. Can be 'up' or 'down'. Default: down
                 
		  For any additional questions, email the author of this script at
		  csanfior@mit.edu'''
         sys.exit()
      elif opt in ("-T"):
         FPKMfile = arg
	 required += [arg]
      elif opt in ("-f"):
         fastaFile = arg
	 required += [arg]
      elif opt in ("-r"):
         cutaway = arg/1.0
      elif opt in ("-R"):
         regType = arg
      elif opt in ("-b"):
         bias = True
      else :
	 print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
         print 'Note: Type \'-h\' for additional help'
	 sys.exit(2)
   args = sys.argv[1:]
   if not args :
	print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if len(required) == 1 :
	print 'Error: Missing argument.' 
	print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(FPKMfile) == False :
	print 'Error:',FPKMfile,'not found' 
	print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(fastaFile) == False :
	print 'Error:',fastaFile,'not found' 
	print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if cutaway.is_integer() == False :
	print 'Error:',cutaway,'is not an interger' 
	print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if regType.lower() not in ['down','up'] :
	print 'Error: -R must be \'down\' or \'up\'.' 
	print 'Usage: TranslationEfficiency.py -T <FPKMTable> -f <fastafile> [optional: -r <retention value> -R <regulation type> ] > MotifsNoRef.fa'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])


########

columns=['Gene Symbol', 'Ribo-Seq: Control (Optimal conditions) FPKM',	'Ribo-Seq: Amino acid starvation FPKM',	'RNA-Seq: Control (Optimal conditions) FPKM', 'RNA-Seq: Amino acid starvation FPKM']

df = pd.read_csv(FPKMfile, sep=" ", skipinitialspace=False, header=None)
df.columns = columns
index = range(len(df)) 
pd.DataFrame(index=index)
newLst = []
transEff = []
for idx in df.index :

	if float(df.ix[:, 'RNA-Seq: Control (Optimal conditions) FPKM'][idx]) - float(df.ix[:, 'RNA-Seq: Amino acid starvation FPKM'][idx]) >= (cutaway/100)*float(df.ix[:, 'RNA-Seq: Control (Optimal conditions) FPKM'][idx]) :
		transEff += ['Transcription Regulation']
		newLst += [0]
	else : 
		result = (float(df.ix[:, 'Ribo-Seq: Control (Optimal conditions) FPKM'][idx]) ) - (float(df.ix[:, 'Ribo-Seq: Amino acid starvation FPKM'][idx]) )
                if float(df.ix[:, 'Ribo-Seq: Control (Optimal conditions) FPKM'][idx]) != 0 :
			transEff += [-result/float(df.ix[:, 'Ribo-Seq: Control (Optimal conditions) FPKM'][idx])]
                else :
                        transEff += [-result] #This fix assumes that the FPKM value for this gene is also small (as to not affect training during SVM implementation)
			
		
		#We will assume that any gene expressed in sufficiently small amounts as to result in 0 FPKM in either control or experimental data sets should be ignored for minimum estimation
		if float(df.ix[:, 'Ribo-Seq: Control (Optimal conditions) FPKM'][idx]) == 0 or float(df.ix[:, 'Ribo-Seq: Amino acid starvation FPKM'][idx]) == 0 :
			newLst += [0]

		else :	
			#To know which gene has the largest difference in translation levels, we will need to divide each substraction by their control FPKM
			newLst += [float(-result)/float(df.ix[:, 'Ribo-Seq: Control (Optimal conditions) FPKM'][idx])]

df['Translational Efficiency'] = transEff

#Minimum or maximum estimation
if regType.lower() == 'down' :
	minID = [i for i,x in enumerate(newLst) if x == min(newLst)]
	minGene = df.ix[:, 'Gene Symbol'][minID].values[0]

elif regType.lower() == 'up' :
	minID = [i for i,x in enumerate(newLst) if x == max(newLst)]
	minGene = df.ix[:, 'Gene Symbol'][minID].values[0]

columns=['Gene Symbol', 'Ribo-Seq: Control (Optimal conditions) FPKM',	'Ribo-Seq: Amino acid starvation FPKM',	'RNA-Seq: Control (Optimal conditions) FPKM', 'RNA-Seq: Amino acid starvation FPKM', 'Translational Efficiency']

df.to_csv('MEP2.fpkm.csv')

#Produce updated csv w/o reference genes
df = df[df['Gene Symbol'] != minGene]
df.reset_index(drop=True,inplace=True)

if bias == True : #Move me up!
	#df = df[df['Ribo-Seq: Control (Optimal conditions) FPKM'] <= 1 ]
	#df = df[df['RNA-Seq: Control (Optimal conditions) FPKM'] <= 1 ] #Fix Me
	#df = df[df['RNA-Seq: Amino acid starvation FPKM'] <= 1 ]
	df.reset_index(drop=True,inplace=True)


 
df.to_csv('MEP2_NoRefGene.fpkm.csv')

#Produce updated fasta

# A nice little file reader
def contentExtractor(files) :
	with open(files) as f:
    		content = f.readlines()
	content = [x.strip('\n') for x in content] 
	return content

content = contentExtractor(fastaFile)
############ Generator 1 (Header and sequence extractor; also removes from fasta gene that will be used as reference)
sequences = []
headers = []
referenceHead = []
referenceSeq = []
for element in content :

	if element[0] == '>' :
		if element != '>'+minGene :
			headers += [element]
			previousElem = ''
		else :
			previousElem = element
			referenceHead += [element]
	else :
		if previousElem != '>'+minGene :
			sequences += [element]
		else :
			referenceSeq += [element]
############

############ Generator 2 (Outputs updated fasta w/o reference genes)
for i in range(len(headers)) :
	print headers[i],'\n',sequences[i]
############

############ Generator 3 (Outputs fasta w/ genes to be used as reference for Needle alignment)
f = open('NeedleReference.fa','w')
for i in range(len(referenceHead)) :
	print >>f, referenceHead[i], '\n', referenceSeq[i]
f.close()
############

#Output MEP2 settings log
f = open('MEP2.settings','w')
print >>f, "[MEP2Settings]"
print >>f, "Regulation Type =", regType.lower()
print >>f, "Retention value =", cutaway
print >>f, "FPKM Table used =",FPKMfile
print >>f, "Fasta used =",fastaFile
print >>f, "Reference Gene =", minGene
f.close()

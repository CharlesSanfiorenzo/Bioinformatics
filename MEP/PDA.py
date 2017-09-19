######Position alignment using matrix. Only use for motif families##############
######If you have any questions, please contact the author of this script#######
######at csanfior@mit.edu ######################################################

import sys, getopt
import pandas as pd
import os.path
import re #Used to retain only numbers from Needle alignment scores


#Note: This module does not take into consideration gaps, as it assumes effective proximity is relative to the final score.

def main(argv):
   global needle
   global matrix
   global fasta
   needle = ''
   matrix = ''
   fasta =''
   required = []
   try:
      opts, args = getopt.getopt(argv,"hN:m:f:")
   except getopt.GetoptError:
      print 'Usage: PDA.py -m <PWM Matrix> -f <sequences(fasta)> > AlignmentScores'
      print 'Note: Type \'-h\' for additional help'
      sys.exit(2)
   for opt, arg in opts:
      if len(sys.argv) == 1 :
	print 'Usage: PDA.py -m <PWM Matrix> -f <sequences(fasta)> > AlignmentScores'
      	print 'Note: Type \'-h\' for additional help'
      	sys.exit(2)
      elif opt == '-h':
         print '''Usage: PDA.py -m <PWM Matrix> -f <sequences(fasta)> > AlignmentScores
                  Author: Charles Sanfiorenzo (2017)
		  --------------------------------------------------------------------
                  Use the Position Dependent Aligner in cases where position in motifs are of upmost importance to
                  their regulatory effect. All motifs must have the same length.Scores will correspond the sequence order in the
		  fasta used.

                  -N | Needle Output
                  -m | Position Weight Matrix. Example of expected input (format of the matrix) for 11 nucleotide-long sequences:
		       0, 1, 0, 0
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 1, 0, 0
                       0, 1, 0, 0
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 0, 0, 1
                       0, 1, 0, 0
                       0, 0, 1, 0
                 
		       Column order represents: A, C, G, T
		       Every row is a position in the sequence (1-11 if the sequence is 11 nucleotides long)
		       Make sure to follow the column order stated above.
	          -f | Sequences (fasta)

		  For any additional questions, email the author of this script at
		  csanfior@mit.edu'''
         sys.exit()
      elif opt in ("-m"):
         matrix = arg
	 required += [arg]
      elif opt in ("-f"):
         fasta = arg
	 required += [arg]
      else :
	 print 'Usage: ApplyPWM.py -N <Needle output> -m <PWM Matrix> -f <sequences(fasta)> > WeightedScores'
         print 'Note: Type \'-h\' for additional help'
	 sys.exit(2)
   args = sys.argv[1:]
   if not args :
	print 'Usage: ApplyPWM.py -N <Needle output> -m <PWM Matrix> -f <sequences(fasta)> > WeightedScores'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if len(required) <= 1 :
	print 'Error: Missing argument.' 
	print 'Usage: ApplyPWM.py -N <Needle output> -m <PWM Matrix> -f <sequences(fasta)> > WeightedScores'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(matrix) == False :
	print 'Error:',matrix,'not found' 
	print 'Usage: ApplyPWM.py -N <Needle output> -m <PWM Matrix> -f <sequences(fasta)> > WeightedScores'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(fasta) == False :
	print 'Error:',fasta,'not found' 
	print 'Usage: ApplyPWM.py -N <Needle output> -m <PWM Matrix> -f <sequences(fasta)> > WeightedScores'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])



#A nice little file reader
def contentExtractor(files) :
	with open(files) as f:
    		content = f.readlines()
	content = [x.strip('\n') for x in content] 
	return content

nonref = contentExtractor(fasta)
sequences = []
headers = []
for element in nonref :
	if element[0] != '>' :
		sequences += [element]
	else :
		headers += [element.strip('>')]

columns = ['A','C','G','T']
df = pd.read_csv(matrix, sep=",", skipinitialspace=False,names = columns)
#print 'Matrix used:\n', df #Prints matrix being used
PositionalPenalty = []
for seq in sequences :
	i = 0
	for base in seq :
		PositionalPenalty += [float(df.ix[:, base][i])]
		i += 1

UpdatedPositionPenalty = [PositionalPenalty[n:n+11] for n in range(0, len(PositionalPenalty), 11)]

#Operation time
weight = []
products = []
for sublst in UpdatedPositionPenalty :
	value = 1
	for element in sublst :
		value += float(element)
	products += [value]



#############

# Applying the weights and outputting...
i = 0
for score in products :
	print 'Score:', float(score)

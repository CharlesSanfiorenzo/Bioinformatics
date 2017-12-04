## OOF, born from a meme. ##################################################################
## OOF finds ALL nested ORFs. ##############################################################
## Any questions about this script, contact the author at csanfior@mit.edu ##
from Bio import SeqIO
from Bio.Alphabet import IUPAC, ProteinAlphabet
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqRecord import SeqRecord
import sys, getopt
import os.path

########Argument parser
def main(argv):
   global src_filename
   global t
   global faa_filename
   src_filename = ''
   t = False
   faa_filename = ''
   required = []
   
   try:
      opts, args = getopt.getopt(argv,"hf:o:t")
   except getopt.GetoptError:
      print 'Usage: OOF.py -f <fasta> -o <outName> [optional: -t]'
      print 'Note: Type \'-h\' for additional help'
      sys.exit(2)
   for opt, arg in opts:
      if len(sys.argv) == 1 :
	print 'Usage: OOF.py -f <fasta> -o <outName> [optional: -t] '
      	print 'Note: Type \'-h\' for additional help'
      	sys.exit(2)
      elif opt == '-h':
         print '''Usage: OOF.py -f <fasta> -o <outName> [optional: -t]
                  Author: Charles Sanfiorenzo (2017)
		  --------------------------------------------------------------------
                  Yes, the name came from a Roblox meme. OOF exhaustively searches for
                  all ORFs found within a given set of sequences.

		  Optional Arguments: 
		  -t | Translates nested ORFs

		  For any additional questions, email the author of this script at
		  csanfior@mit.edu'''
         sys.exit()
      elif opt in ("-f"):
         src_filename = arg
	 required += [arg]
      elif opt in ("-o"):
         faa_filename = arg
	 required += [arg]
      elif opt in ("-t"):
         t = True

      else :
	 print 'Usage: OOF.py -f <fasta> -o <outName> [optional: -t]'
         print 'Note: Type \'-h\' for additional help'
	 sys.exit(2)
   args = sys.argv[1:]
   if not args :
	print 'Usage: OOF.py -f <fasta> -o <outName> [optional: -t]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if len(required) == 1 :
	print 'Error: Missing argument.' 
	print 'Usage: OOF.py -f <fasta> -o <outName> [optional: -t]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(src_filename) == False :
	print 'Error:',src_filename,'not found' 
	print 'Usage: OOF.py -f <fasta> -o <outName> [optional: -t]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()


if __name__ == "__main__":
	main(sys.argv[1:])

countNope=0
codonTableID = 1
n = 3
stopPosList = []
output_handle = open(faa_filename, "w")
output_handle2 = open(str(faa_filename)+"Longest", "w")
input_handle  = open(src_filename, "r")

if t == False :
	def find_best_frame(seq, directionsToConsider="forward", tranlationTable=1):

	    allPossibilities = []
	    if directionsToConsider in ("forward"):
	        # start translation from 1, 2 and 3 nucleotide
	        for frame in range(len(seq)):
		    #if it starts with M (add this)
		    if str(seq[frame:][:3]) == 'ATG' :
			frame_split = [seq[frame:].strip()[i:i+n] for i in range(0,len(seq[frame:].strip()),n)]
			#print Seq(frame_split)
			if "TAA" in frame_split or "TGA" in frame_split or "TAG" in frame_split : #if there is a stop codon only index until stop codon
				#stopPosList = [i for i,x in enumerate(frame_split) if x in ["TGA","TAG","TAA"]] 
				stopPosList = []
				count = 0
				for i in frame_split :
					if i in ["TGA","TAG","TAA"] :
						stopPosList += [count]
					count += 1
					#print i
					
				stopPos = (stopPosList[0]*3+3) #Get true stop position
				no_trans = str(seq[frame:frame+stopPos])
				#CHECKERS: Comment them
				#print stopPosList
				#print stopPos
				#print seq[frame:]
				#print seq[frame:frame+15]
			else : 
				#TRUNCATION STEP
				if len(seq[frame:]) % 3 == 0 :
	            			no_trans = str(seq[frame:])
				else :
					#print len(frame_split[-1]) 
					no_trans = str(seq[frame:len(seq)-len(frame_split[-1])])
	            	allPossibilities.append(no_trans)
	
	    # Choose the longest peptide from an ORF that starts with an ATG codon
	    lenPep = [len(i) for i in allPossibilities]
	    bestFrame = allPossibilities[lenPep.index(max(lenPep))]
	    #If length of protein > 10 a.a
	    #if len(bestFrame) >= 10*3 :
	    return Seq(bestFrame, alphabet=IUPAC.ambiguous_dna)

	def find_all_frame(seq, directionsToConsider="forward", tranlationTable=1):

	    allPossibilities = []
	    if directionsToConsider in ("forward"):
	        # start translation from 1, 2 and 3 nucleotide
	        for frame in range(len(seq)):
		    #if it starts with M (add this)
		    if str(seq[frame:][:3]) == 'ATG' :
			frame_split = [seq[frame:].strip()[i:i+n] for i in range(0,len(seq[frame:].strip()),n)]
			#print Seq(frame_split)
			if "TAA" in frame_split or "TGA" in frame_split or "TAG" in frame_split : #if there is a stop codon only index until stop codon
				#stopPosList = [i for i,x in enumerate(frame_split) if x in ["TGA","TAG","TAA"]] 
				stopPosList = []
				count = 0
				for i in frame_split :
					if i in ["TGA","TAG","TAA"] :
						stopPosList += [count]
					count += 1
					#print i
					
				stopPos = (stopPosList[0]*3+3) #Get true stop position
				no_trans = str(seq[frame:frame+stopPos])
				#CHECKERS: Comment them
				#print stopPosList
				#print stopPos
				#print seq[frame:]
				#print seq[frame:frame+15]
			else : 
				#TRUNCATION STEP
				if len(seq[frame:]) % 3 == 0 :
	            			no_trans = str(seq[frame:])
				else :
					#print len(frame_split[-1]) 
					no_trans = str(seq[frame:len(seq)-len(frame_split[-1])])
	            	allPossibilities.append(no_trans)
	
	    # Choose the longest peptide from an ORF that starts with an ATG codon
	    lenPep = [len(i) for i in allPossibilities]
	    bestFrame = allPossibilities[lenPep.index(max(lenPep))]
	    #If length of protein > 10 a.a
	    #if len(bestFrame) >= 10*3 :
	    return Seq(allPossibilities, alphabet=IUPAC.ambiguous_dna)
	
	for seq_record in SeqIO.parse(input_handle, "fasta", alphabet=IUPAC.ambiguous_dna):
	    try:
	        bestFrame = find_best_frame(seq_record.seq,"forward", codonTableID)
	        bestFrameRecord = SeqRecord(bestFrame, seq_record.name)
	        bestFrameRecord.description = seq_record.description
	
	        SeqIO.write(bestFrameRecord, output_handle2, "fasta")

##########################################################################################
		bestFrame = find_all_frame(seq_record.seq,"forward", codonTableID)
		print bestFrame
	        bestFrameRecord = SeqRecord(bestFrame, seq_record.name)
	        bestFrameRecord.description = seq_record.description
	
	        SeqIO.write(bestFrameRecord, output_handle, "fasta")

	    except Exception as inst:
		print "No ORF found here. Moving on"

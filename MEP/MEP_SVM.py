######################MEP Support Vector Machine.##############################
######If you have any questions, please contact the author of this script#######
######at csanfior@mit.edu ######################################################

import numpy as np
import sys, getopt
import pandas as pd
import os.path
import random
from sklearn import svm, linear_model
from sklearn.model_selection import KFold, cross_val_score, GridSearchCV #Change 'model_selection' to 'cross_validation' for older versions of sklearn
import matplotlib.pyplot as plt
import ConfigParser

########Argument parser
def main(argv):
   global FPKMfile
   global db
   global n
   global v
   global randomizer
   global over
   global threads
   threads = 1
   FPKMfile = ''
   db = ''
   n = 40.0
   v = ''
   over = False
   randomizer = False
   required = []
   try:
      opts, args = getopt.getopt(argv,"hd:M:v:n:@:",["randomizer","over-fit"])
   except getopt.GetoptError:
      print 'Usage: MEP.py -d <MEP database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
      print 'Note: Type \'-h\' for additional help'
      sys.exit(2)
   for opt, arg in opts:
      if len(sys.argv) == 1 :
	print 'Usage: MEP.py -d <MEP database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
      	print 'Note: Type \'-h\' for additional help'
      	sys.exit(2)
      elif opt == '-h':
         print '''Usage: MEP.py -d <MEP database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]
		  Example: TranslationalEfficiency.py -d MEP2.db -M MEP2_NoRefGene.fpkm.csv -n 40
                  Author: Charles Sanfiorenzo (2017)
		  --------------------------------------------------------------------
                  Running this script will run MEP's Support Vector & Linear Regression models; the most
                  accurate model for each will be outputted and graphed, along with a table containing 
                  accuracy comparison between models.

                  -d | MEP2 database produced by DBCreator.py
                  -M | MEP2 FPKM Table w/o Reference Gene. It is used as a 'sanity' check. Produced by TranslationEfficiency.py
                  -v | Validation (testing) data set (copy the format of MEP.db, but add 'Label' as a column; see examples on GitHub). 
                       If validation genes are within MEP.db, supply this option with a list of corresponding Gene Symbols and an assigned label 
                       (1 for effect observed, -1 for no effect observed).
                  
                  Optional arguments:
                  -n | number of genes to be used for training when running MEP. Default: 40
                       Example: for -n 40, 40 genes will be chosen for positive labels,
                       and 40 genes will be chosen for negative labels (a total of 80 genes 
                       will be used to train).
	--randomizer | selection of test sample is random (not recommended!) Default testing: cross-validation
          --over-fit | Purposely overfit the SVM model by using Translation Efficiency as a feature.
                 
		  For any additional questions, email the author of this script at
		  csanfior@mit.edu'''
         sys.exit()
      elif opt in ("-d"):
         db = arg
	 required += [arg]
      elif opt in ("-M"):
         FPKMfile = arg
	 required += [arg]
      elif opt in ("-v"):
         v = arg
         required += [arg]
      elif opt in ("-n"):
         n = arg/1.0
      elif opt in ("--randomizer"):
         randomizer = True
      elif opt in ("--over-fit"):
         over = True
	 randomizer = True
      elif opt in ("-@"):
         threads = arg
      else :
	 print 'Usage: MEP2.py -d <MEP2 database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
         print 'Note: Type \'-h\' for additional help'
	 sys.exit(2)
   args = sys.argv[1:]
   if not args :
	print 'Usage: MEP2.py -d <MEP2 database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if len(required) <= 2 :
	print 'Error: Missing argument.' 
	print 'Usage: MEP2.py -d <MEP2 database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(FPKMfile) == False :
	print 'Error:',FPKMfile,'not found' 
	print 'Usage: MEP2.py -d <MEP2 database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(db) == False :
	print 'Error:',db,'not found' 
	print 'Usage: MEP2.py -d <MEP2 database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if os.path.isfile(v) == False :
	print 'Error:',v,'not found' 
	print 'Usage: MEP2.py -d <MEP2 database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if n.is_integer() == False :
	print 'Error:',n,'is not an interger' 
	print 'Usage: MEP2.py -d <MEP2 database> -M <MEP_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()
   if threads.is_integer() == False :
	print 'Error: Number of threads must be an interger.' 
	print 'Usage: MEP2.py -d <MEP2 database> -M <MEP2_NoRefGene.fpkm.csv> -v <validation set> [optional: -n <label number>]'
        print 'Note: Type \'-h\' for additional help'
	sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])


########
config = ConfigParser.RawConfigParser(allow_no_value=True)
config.read("MEP2.settings")
#config.get()


#A nice little file reader
def contentExtractor(files) :
	with open(files) as f:
    		content = f.readlines()
	content = [x.strip('\n') for x in content] 
	return content

#Load settings
settings = contentExtractor('MEP2.settings')

#Load MEP2 database
df = pd.read_csv(db, header=0, sep=",", index_col=0)


#Load MEP2 FPKM Table
df2 = pd.read_csv(FPKMfile, header=0, sep=",", index_col=0)
columns=['Gene Symbol', 'Ribo-Seq: Control (Optimal conditions) FPKM', 'Ribo-Seq: Amino acid starvation FPKM', 'RNA-Seq: Control (Optimal conditions) FPKM', 'RNA-Seq: Amino acid starvation FPKM', 'Translational Efficiency']
df2.columns = columns
df2 = df2[df2['Translational Efficiency'] != 'Transcription Regulation'] #Eliminates genes that showed Transcription Regulation
df = df[df['Translational Efficiency'] != 'Transcription Regulation'] #Eliminates genes that showed Transcription Regulation
df.reset_index(drop=True,inplace=True) #Reset index for df
df2.reset_index(drop=True,inplace=True) #Reset index for df2


######## Generator 1 (a few calculations for the inclusion of outliers)
newLst = []
for idx in df2.index :
	result = (float(df2.ix[:, 'Ribo-Seq: Control (Optimal conditions) FPKM'][idx]) ) - (float(df2.ix[:, 'Ribo-Seq: Amino acid starvation FPKM'][idx]) )

	#We will assume that any gene expressed in sufficiently small amounts as to result in 0 FPKM in either control or experimental data sets should be ignored for minimum estimation
	if float(df2.ix[:, 'Ribo-Seq: Control (Optimal conditions) FPKM'][idx]) == 0 or float(df2.ix[:, 'Ribo-Seq: Amino acid starvation FPKM'][idx]) == 0 :
		if settings[1] == 'Regulation Type = down' :
			newLst += [999999] #outlier
		elif settings[1] == 'Regulation Type = up' : 
			newLst += [-999999] #outlier
	else :	
		#To know which gene has the largest difference in translation levels, we will need to divide each substraction by their control FPKM
		newLst += [float(-result)/float(df2.ix[:, 'Ribo-Seq: Control (Optimal conditions) FPKM'][idx])]

#######

trainDF = pd.DataFrame() #Training data set
columns=['Gene Symbol', 'Sequence', 'Alignment Score', 'Translational Efficiency', 'label']
trainDF.columns
trainDF['label'] = np.nan

###### Minimum or maximum estimarion estimation (label)
newLst = [float(x) for x in newLst]
if settings[1] == 'Regulation Type = down' :
	newLstSrt = sorted(newLst,key=float)[0:int(n)]

elif settings[1] == 'Regulation Type = up' :
	newLstSrt = sorted(newLst,key=float)[-1:int(-n)]

newLstRnd = []
for element in newLstSrt :
	newLstRnd += [round(element,4)]

df['Translational Efficiency'] = df['Translational Efficiency'].astype('float64') #Computer Science magic
df = df.round(decimals=4) #Computer Science magic

if over == True : #This one eliminates from df, condition it
	for element in newLstRnd :
		trainDF = trainDF.append(df[df['Translational Efficiency'] == element])
		trainDF.fillna(value=1.0, method=None, inplace=True)
		df = df[df['Translational Efficiency'] != element] #Eliminates genes that showed Transcription Regulation
		df.reset_index(drop=True,inplace=True) #Reset index for df
	trainDF.reset_index(drop=True,inplace=True) #Reset index for trainDF
	######
                
	###### Closest to 0 (label)
	newLstSrt = [i for d, i in sorted((abs(x-0), x) for x in newLst)[:int(n)]] # Get n closest numbers to 0
	newLstRnd = []
	for element in newLstSrt : #Which should only contain 40 genes by default
		newLstRnd += [round(element,4)]

	for element in newLstRnd :
		trainDF = trainDF.append(df[df['Translational Efficiency'] == element])
		trainDF.fillna(-1.0, inplace=True)
		df = df[df['Translational Efficiency'] != element] #Eliminates genes that showed Transcription Regulation
		df.reset_index(drop=True,inplace=True) #Reset index for df
	trainDF.reset_index(drop=True,inplace=True) #Reset index for trainDF

	######

###### Test Data
#Randomizer for 80 (not recommended!)
if over == True :
	testLst = []
	for idx in range(80) :
		testLst += [random.randint(0, len(df.index))]
	testDF = df.ix[testLst]
	testDF.reset_index(drop=True,inplace=True)

	X_test = [] #test features
	y_test = [] #test labels

	for idx in testDF.index :
		X_test += [[testDF.ix[:, 'Alignment Score'][idx]]]
		X_test[idx] += [(float(testDF.ix[:,'Translational Efficiency'][idx]))]

		if (float(testDF.ix[:,'Translational Efficiency'][idx])) < 0 :
			y_test += [1]
		else :
			y_test += [-1]

########SVM Implementation
if over == True :
	#Training feature & label extraction
	X_train = [] # features
	y_train = [] # labels
	for idx in trainDF.index :
		X_train += [[trainDF.ix[:, 'Alignment Score'][idx]]]
		X_train[idx] += [(float(trainDF.ix[:,'Translational Efficiency'][idx]))]
		y_train += [int(trainDF.ix[:, 'label'][idx])]

elif over == False : #If number of one label is several times larger than the other, use samplified stratification #REMEMBER TO CHANGE THE LABEL
	content = contentExtractor('MEP2.settings')
	print df #REMOVE
	X = []
	y = []
	for idx in df.index :
		X += [[df.ix[:, 'Alignment Score'][idx]]]
		X[idx] += [(float(df.ix[:,str(content[7][10:])][idx]))]
		if (float(df.ix[:,'Translational Efficiency'][idx])) < 0 : #CHANGE ME
			y += [1]
		else :
			y += [-1]
	#Kfold crossvalidation (Default)
	kf = KFold(random_state=None, shuffle=False)
	for train_index, test_index in kf.split(X) :
		X_train, X_test = np.array(X)[train_index], np.array(X)[test_index]
		y_train, y_test = np.array(y)[train_index], np.array(y)[test_index]
	
#Kernel & parameter selector
kernels=['linear','rbf','sigmoid']
kernelLog = []
GeneralAccuracyLst = []
bestMEP2 = []
for kernel in kernels :
	MEP2baby = svm.SVC(kernel=kernel)
	c_range = range(1,4,10)
	MEP2 = GridSearchCV(MEP2baby, param_grid=dict(C=c_range)) #Exhaustive grid search
	MEP2.fit(X_train, y_train)
	kernelLog += [kernel]

		#Crossvalidation (Default)
	if over == False :
		GeneralAccuracy = cross_val_score(MEP2, X, y, cv=kf, n_jobs = 2)
		print GeneralAccuracy #REMOVE
		#Random sampling for validation
	elif over == True:
		GeneralAccuracy = MEP2.score(X_test, y_test)*100  #Validation step
			
	bestMEP2 += [MEP2]
	GeneralAccuracyLst += [GeneralAccuracy]

if over == False :

	Means = []	
	for i in range(len(GeneralAccuracyLst)) :
		Means += [GeneralAccuracyLst[i].mean()]
	maxID = [i for i,x in enumerate(Means) if x == max(Means)][0]
	print "[MEP2 SVM Results]"
	print "Best Kernel:", kernelLog[maxID]
	print("Accuracy values:")
	i=0
	for accuracies in GeneralAccuracyLst :
		print accuracies, kernelLog[i] 
		i += 1
	print("Mean Accuracy: %0.2f (+/- %0.2f)" % (max(Means), max(Means).std() * 2)) #Fix std
	print "Best Accuracy:", max(GeneralAccuracyLst[maxID])

elif over == True :
	maxID = [i for i,x in enumerate(GeneralAccuracyLst) if x == max(GeneralAccuracyLst)][0]
	print "Best Kernel:", kernelLog[maxID]
	print "General Accuracy Score: ", max(GeneralAccuracy)  #Validation step


#Once the best kernel is chosen...
#MEP2.predict([[2., 2.]]) #Prediction step

###Plotting time baby

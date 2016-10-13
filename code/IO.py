
'''
Created on 2010
@author: Panagiotis Achlioptas
@contact: pachlioptas@gmail.com
@copyright: You are free to use, change, or redistribute this code in any way you want 
for non-commercial purposes. 
'''

import sys, cPickle, numpy
    
def loadSNPs(genoFile, verbose = True):
    snps = []
    if verbose:
        print("Loading SNPs from file \"" + genoFile +"\"")
    with open(genoFile) as snpFile:
        for counter, content in enumerate(snpFile):
            snp = [int(i) for i in content.split(None)]                 

            if len(set(snp)) == 1:
                if verbose:
                    print "SNP " + str(counter) + " has zero variance, so we discard it."
                
            else:
                snps.append(snp)

    numOfSnps   = len(snps)
    individuals = len(snps[0])
    if verbose:
        print "\nThe number of SNPs loaded is " + str(numOfSnps) + "."
        print "Each SNP has values for", individuals, "individuals.\n"
    
    
    for i in xrange(numOfSnps):
        for j in xrange(individuals):            
            assert (snps[i][j] == 0 or snps[i][j] == 1 or snps[i][j] == 2) #Sanity Check
        
    snps = numpy.array(snps, dtype=numpy.int8)

    return snps



def promptUserForInput(verbose = True):

    if len(sys.argv) != 8 :
        print ('Incorrect Input Arguments Were Given. Usage:\n\n'\
               '1st   Cases Genotypes\n'\
               '2nd   Controls Genotypes\n'\
               '3rd   Lower  Adaptive Bound\n'\
               '4rt   Upper  Adaptive Bound\n'\
               '5th   Number Of Repetitions\n'\
               '6th   Number of Hyperplanes\n'\
               '7th   Number of Top SNPs to report\n'
               )
    
        print sys.argv
        sys.exit()

    try:
        lowerAdaptB         = int(sys.argv[3])
        upperAdaptB         = int(sys.argv[4])
        repetitions         = int(sys.argv[5])
        hyperplanesNum      = int(sys.argv[6])
        shouldReport        = int(sys.argv[7])

    except Exception as inst:
        print type(inst)
        sys.exit()
    
    cases = loadSNPs(sys.argv[1], verbose)
    ctrls = loadSNPs(sys.argv[2], verbose)
    
    if verbose:
        print "--Input PARAMETERS--"
        print "Lower adaptive bound " + str(lowerAdaptB) + "."
        print "Upper adaptive bound " + str(upperAdaptB) + "." 
        print "Number of repetitions " + str(repetitions) + "."
        print "Number of Hyperplanes " + str(hyperplanesNum) + "."
        print "Number of top-Pairs that will be kept " + str(shouldReport) + "."
        
    return cases, ctrls, lowerAdaptB, upperAdaptB, repetitions, hyperplanesNum, shouldReport
    

def saveData(fileName, *args):    
    myFile = open(fileName, "w")    
    cPickle.dump(len(args), myFile)
    for item in args:
        cPickle.dump(item, myFile)
    myFile.close()    



def loadSavedData(fileName):
    inFile = open(fileName, "r")
    size = cPickle.load(inFile)
    for i in xrange(size):
        yield cPickle.load(inFile)        
    inFile.close()
    


def numericStringToIntegerList(inString, outLength):
    res = []
    for i in xrange(outLength):
        res.append(int(inString[i]))
    return res


def storeSNPsAsText (file, snps, metadata = None):
    '''
    Preconditions: len(metadata) == len(snps), aka each SNP carries one set of meta-data for it
    Each SNP with its associated meta data will be store in a single line. Values are space separated.
    '''
    if metadata != None:
        with open (file, "w") as fout:
            for meta, snp in zip(metadata, snps):
                for m in meta:
                    fout.write(m)
                    fout.write(" ")
                for s in snp:
                    fout.write(str(s))
                    fout.write(" ")
                fout.write("\n")
    else:
        with open (file, "w") as fout:
            for snp in snps:                
                for s in snp:
                    fout.write(str(s))
                    fout.write(" ")
                fout.write("\n")



##
## LOADING METHODS FOR PARTICULAR INPUT SNP FILES
##
def load_HAP_SNPs(genotypeFile, numberOfSNPs, individuals, selectionMethod = "continuous", randomGroup = None, verbose = False):

    snps = []
    metadata = []

    with open(genotypeFile) as snpsFile:    
        
        if selectionMethod == "continuous":            
            i = 0
            for line in snpsFile:                
                if i < numberOfSNPs:
                    snps.append(line.split(None))
                i += 1                
            
        elif selectionMethod == "random":
            allLines = snpsFile.readlines()
            print "The Input Genotype File has %d SNPs" % len(allLines)
            i = 0            
            if randomGroup == None:
                import random
                randomGroup = random.sample(xrange(len(allLines)), numberOfSNPs)
                              
            for line in allLines:
                if i in randomGroup:
                    snps.append(line.split(None))
                i += 1
    
        else: assert(False)
    
    assert(len(snps) == numberOfSNPs)
    
    mtd = 4  #they keep 4 meta data values for each SNP 
    
    for snp in snps:        
        metadata.append(snp[:mtd])        
        snp[:] = snp[mtd:individuals+mtd] 
        snp[:] = numericStringToIntegerList(snp, individuals)        

        assert(len(snp) == individuals)
        for value in snp:
            assert(type(value) == int and (value == 0 or value == 1 or value == 2))
                
    return snps, metadata, randomGroup 

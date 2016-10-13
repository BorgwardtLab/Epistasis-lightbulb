'''
Created on 2010
@author: Panagiotis Achlioptas
@contact: pachlioptas@gmail.com
@copyright: You are free to use, change, or redistribute this code in any way you want 
for non-commercial purposes. 
'''

import numpy, copy, warnings, operator, scipy.stats as stats


def flipVectorInArray(myArray, linePos):
    columns = len(myArray[0])
    for i in xrange(columns):
        if myArray[linePos][i] == 0: myArray[linePos][i] = 1
        else                       : myArray[linePos][i] = 0
        

def flipArray(myArray):
    '''
    Returns a flipped version of the binary input matrix "MyArray".
    The input matrix remains unchanged.
    '''
    arrayCopy = copy.deepcopy(myArray) 
    for i in xrange(len(myArray)):
            flipVectorInArray(arrayCopy, i)
    
    return arrayCopy


         
def insertInSortedList(item, myList):  
    '''
    The value with witch the sorting is done is the first entry of each item in the input -myList-, i.e. myList[0][i].
    The -myList- must be sorted in a descending fashion, i.e. myList[0][i] <= myList[0][j] when i<j
    '''
    if item[0] < myList[0][0] : return    
    
    myList[0] = item        
    for i in xrange(len(myList) - 1):
        if myList[i][0] > myList[i+1][0]: #swap
            temp = myList[i+1]
            myList[i+1] = myList[i]
            myList[i] = temp  
        else: return


def centerMatrixRows( myMatrix ):
    centered = []
    for row in myMatrix:
        meanValue = numpy.mean(row)        
        centered.append([i-meanValue for i in row])
    return centered
#        for i in xrange(len(row)):  why this doesn't work?
#            row[i] = row[i] - meanValue



def sideOfHyperPlane(a):
    a[a >= 0] = 1
    a[a <  0] = 0
    a = a.astype('int8')
    return a



def hypeplanesHashing(sub1, sub2, hyperplanesNum):
    '''
        Generate a collection of hyper-planes to hash the input matrices -sub1, sub2-.
        Sub1 and Sub2 must have the same number of rows (SNPs) but may vary at the number of their columns (individuals).
        If such an imbalance exists we start by hashing the smallest of the two matrices and we expand the collection
        of hyper-planes we used, by adding to it as many dimensions as needed to couple/hash the larger input matrix.
        
        -HyperplanesNum- controls the number of hyper-planes used, thus also the length of the output SNP (row) vectors.
    '''
    
    sub1Indivs = len(sub1[0])
    sub2Indivs = len(sub2[0])
    
    diff       = abs(sub1Indivs - sub2Indivs)
    
    sub1 = centerMatrixRows(sub1)
    sub2 = centerMatrixRows(sub2)
    
    for row1, row2 in zip(sub1, sub2): 
        assert(numpy.mean(row1) < 0.001) #due to numerical issues, the centering will not be exact.   
        assert(numpy.mean(row2) < 0.001) #we assert that it is at least 'close' to zero.

    
    if sub1Indivs < sub2Indivs:
        hyperplanesSmall     = numpy.random.randn(sub1Indivs, hyperplanesNum)
        sub1Hashed           = numpy.dot(sub1, hyperplanesSmall)
        hyperplanesDif       = numpy.random.randn(diff, hyperplanesNum)
        hyperplanesLarge     = numpy.vstack((hyperplanesSmall, hyperplanesDif))
        sub2Hashed           = numpy.dot(sub2, hyperplanesLarge)

    elif sub2Indivs < sub1Indivs:
        hyperplanesSmall     = numpy.random.randn(sub2Indivs, hyperplanesNum)
        sub2Hashed           = numpy.dot(sub2, hyperplanesSmall)
        hyperplanesDif       = numpy.random.randn(diff, hyperplanesNum)
        hyperplanesLarge     = numpy.vstack((hyperplanesSmall, hyperplanesDif))
        sub1Hashed           = numpy.dot(sub1, hyperplanesLarge)

    else :
        hyperplanes          = numpy.random.randn(sub1Indivs, hyperplanesNum)        
        sub1Hashed           = numpy.dot(sub1, hyperplanes)
        sub2Hashed           = numpy.dot(sub2, hyperplanes)

    assert  sub1Hashed.shape == sub2Hashed.shape
        
    sub1            = sideOfHyperPlane(sub1Hashed) #Convert to Binary: encode the side of the hyperplane
    sub2            = sideOfHyperPlane(sub2Hashed) #
            
    return sub1, sub2




def proposePairs(cases, ctrls, snpsHashed, snpsHashedPerturbed, sampleSize, bucketSizeBounds, shouldReport, maxTimes):   
    
    if shouldReport == 1:
        warnings.warn("Using \"insertInSortedList()\" is very inefficient when you are looking for a single causal SNP pair.")     
    
    numOfSnps             = len(snpsHashed)
    individuals           = len(snpsHashed[0])                      # Number of dimensions for cases||ctrls after hashing.
    sparse                = bucketSizeBounds[0]                     # Bounds for adaptive sampling.
    bloated               = bucketSizeBounds[1]                     #
    
    pairsEncountered      = set()                                   # Track pairs of SNPs to measure correlation only once.
    topKpairs             = [(-1,-1) for i in xrange(shouldReport)] # Best pairs found. 
    
    assert(numOfSnps  == len(snpsHashedPerturbed) == len(cases) == len(ctrls))
        
    
    sample_matrix1        = numpy.empty((individuals, numOfSnps), dtype=numpy.int8)    
    sample_matrix2        = numpy.empty((individuals, numOfSnps), dtype=numpy.int8)
    allVotes              = 0
    convergenceOfSampling = []                                      # Keep sample size per iteration (for post-statistics).
     
    for i in xrange(maxTimes):
        pairsInRound = 0
        buckets      = dict()    
                
        for i in xrange(sampleSize):
            random_individual    = numpy.random.randint(individuals)                        
            sample_vector        = snpsHashed[:, random_individual]
            sample_matrix1[i, :] = sample_vector
            sample_vector        = snpsHashedPerturbed[:, random_individual]
            sample_matrix2[i, :] = sample_vector
                            
        for i in xrange(numOfSnps):
            current_snp = tuple(sample_matrix1[0:sampleSize, i])        
            if current_snp not in buckets:
                buckets[current_snp] = [i]                    
            else:                        
                buckets[current_snp].append(i)
                        
        for i in xrange(numOfSnps):
            current_snp = tuple(sample_matrix2[0:sampleSize, i])        
            if current_snp in buckets:
                for j in buckets[current_snp]:            
                    if i < j:
                        pairsInRound += 1
                        if (i,j) not in pairsEncountered:                                                                                                      
                            r1 = stats.pearsonr(cases[i], cases[j])[0]
                            r2 = stats.pearsonr(ctrls[i], ctrls[j])[0]                        
                            ro = abs(r1 - r2)                                                
                            pairsEncountered.add((i,j))                        
                            if ro > topKpairs[0][0]: insertInSortedList((ro, (i,j)), topKpairs)                     
                                                                                                 
        
        if pairsInRound > bloated:
            sampleSize +=1
        elif pairsInRound < sparse:
            sampleSize -=1
        assert(sampleSize > 1)
        convergenceOfSampling.append(sampleSize)
        allVotes += pairsInRound
        
    avSamplingSize = int(sum(convergenceOfSampling)/float(len(convergenceOfSampling)))    
    topKpairs.reverse()
    
    diffPairsEnc = len(pairsEncountered) 
    del pairsEncountered
    return topKpairs, allVotes, avSamplingSize, diffPairsEnc





def evaluateOurSolution(BFsolution, ourSolution, typeOfSol = "no_counters"):
    '''
    Given a solution (a set of pairs) that our algorithm discovered and the ground truth, returns an evaluation of our solution.
    
    BFsolution: A list with the set of pairs that constitute the optimal solution. The list is in descending order of importance, e.g.
    correlation_diference of BFsolution(k) > BFsolution(k+1).
    
    Returns: a list comprised by each pair of the -ourSolution- that is also a pair in the BFsolution. For such a pair we return
    also its position in the BFsolution in a tuple: (pair, position).
    '''    
    
    if typeOfSol == "classic":
        index = 0
    elif typeOfSol == "no_counters":
        index = 1
    else:
        assert(False)
    
    bestPairsBfDict = dict()
    
    for i, pair in enumerate(BFsolution):
        bestPairsBfDict[pair[0]] = i
                                
    foundings = list()
    if typeOfSol == "no_counters":          
        previousRank = -1                   
        previousCorr = ourSolution[0][0]

    for i, sol in enumerate(ourSolution):
        
        pair  = sol[index]        #correlation of (x,y) pair must be equal to this of (y,x)
        pairT = (sol[index][1], sol[index][0])
        
        if pair in bestPairsBfDict:
            
            assert(pairT not in bestPairsBfDict)                        
            rank = bestPairsBfDict[pair]
            foundings.append((i, rank))
                                    
            if typeOfSol == "no_counters": #sanity check
                assert (previousRank < rank or previousCorr == sol[0])
                previousRank = rank
                previousCorr = sol[0]

        if pairT in bestPairsBfDict:
            assert(pair not in bestPairsBfDict)                        
            rank = bestPairsBfDict[pairT]
            foundings.append((i, rank))                        
                        
            if typeOfSol == "no_counters": #sanity check
                assert (previousRank < rank or previousCorr == sol[0])
                previousRank = rank
                previousCorr = sol[0]
    
    foundings.sort(key = operator.itemgetter(1))    # This is necessary because we have pairs that have the same correlation
                                                    # and their relative position in the BFsolution is random.
    return foundings





def recall_combined(solutions1, solutions2, topK):    
    combined = set()  
    
    stop1 = min(len(solutions1), topK) 
    stop2 = min(len(solutions2), topK)
        
    for i in xrange(stop1):    
#        if solutions1[i][0] > lookingRange : break   #after sorting the values in evalSolution don't need this.
        if solutions1[i][1] < topK : combined.add(solutions1[i][1] ) 
        
    for i in xrange(stop2):    
#        if solutions2[i][0] > lookingRange : break
        if solutions2[i][1] < topK : combined.add(solutions2[i][1] )

    return len(combined)




def kPairsWithMaxDiffPearCorr(k, matrix1, matrix2):
    '''
        Find in a Brute Force fashion the -k- pairs formed by 2 rows of the input matrices, that have the maximum difference
        in Pearson's correlation, between -matrix1- and -matrix2- .
    '''

    assert(len(matrix1)== len(matrix2))    

    bestPairs = [(-1,-1) for i in xrange(k)]
    
    for i in xrange(len(matrix1)):    
        for j in xrange(i+1, len(matrix1)):                        
                                
            cor1 = stats.pearsonr(matrix1[i], matrix1[j])[0]
            cor2 = stats.pearsonr(matrix2[i], matrix2[j])[0]                        
            diff = abs(cor1 - cor2)
 
#            item = (diff, (i,j), cor1, cor2)    #you could keep also the correlation of each pair
            if diff > bestPairs[0][0]: 
                item = (diff, (i,j))
                insertInSortedList(item, bestPairs)            

    #reverse the list in order the sorting to be descending
    revBestPairs = list()
    for i in xrange(1,k):
        revBestPairs.append((bestPairs[-i][1], bestPairs[-i][0]))
    revBestPairs.append((bestPairs[0][1], bestPairs[0][0]))
    
    return revBestPairs
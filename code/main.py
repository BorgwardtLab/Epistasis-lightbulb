'''
Created on 2010
@author: Panagiotis Achlioptas
@contact: pachlioptas@gmail.com
@copyright: You are free to use, change, or redistribute this code in any way you want 
for non-commercial purposes. 


Quick review of the steps we follow to find the causal SNPs based on the algorithm that we proposed in the paper:
            "Two-Locus Association Mapping in Subquadratic Time, P. Achlioptas et al., KDD 2011" 

Input: a set of -n- ternary valued SNPs associated with:
        -m1- individuals (a matrix of "Cases" with size n X m1) 
  and   -m2- individuals (a matrix of "Controls" with size n X m2).

Goal: find the -k- pairs of SNPs for which the difference of their Pearson Correlation as measured in the Cases and the Controls matrices
is maximized. 
 
The steps of the subquadratic approximation algorithm that we proposed are:

1. Center all the SNPs around zero (force them to have zero mean). The SNPs are now Real-valued vectors. 
2. Form a random collection of -lambda- hyperplanes (each coordinate is a Standard Gaussian r.v.).
3. Convert each SNP (in each matrix) to a -lambda- dimensional binary vector by determining its sign with respect to each hyperplane.
   Now we have two new binary matrices of the same size: n X lambda, formed by the hashed SNPs. Call theses matrices Cases* and Controls* respectively.
4a. Form a new matrix M1 of size: n X 2*lambda, by column-concatenating Cases* with the complement of Controls*.
4b. Form a new matrix M2 of size: n X 2*lambda, by column-concatenating the complement of Cases* with Controls*.
5. Apply the algorithm ProposePairs to the matrix M1 and then to M2.
6. Merge the results from step 5.

Algorithm ProposePairs:
    Input: 'snpMatrix', 'lowerAdaptiveBound', 'higherAdaptiveBound', 'shouldReport', 'repetitions'
        
        -bestPairsFound- = A Heap of size 'shouldReport' (keeps the discovered SNP pairs with the maximum correlation found)
        -s-              = 1/4 of the number of columns of 'snpMatrix'  (initial size of samples taken from each SNP) 
        
    Repeat 'repetitions' times:
        Form a new submatrix -snpMatrixSampled- by sampling -s- columns of the 'snpMatrix'.
        
        Hash all the rows of the -snpMatrixSampled-. If a pair of rows end up in the same bucket (i.e. the pair is comprised by two identical vectors),
        then, measure the Pearson Correlation of the corresponding SNPs (this pair is a "proposed" pair).
        
        Insert the "proposed" pair in the BestPairsFound if its correlation is larger than the smallest residing there.
        
        If the number of "proposed" pairs is less than 'lowerAdaptiveBound'
            s -= 1
        If the number of "proposed" pairs is more than 'upperAdaptiveBound'
            s += 1
'''



import time, sys, IO
from methods import *


if __name__ == '__main__':
    try: 
        import psyco
        psyco.full()
    except:
        print "\n Psyco is a module that speeds up Python's code, it's usage is optional. Psyco is not available.\n"


    cases, ctrls, lowerAdaptB, upperAdaptB, repetitions, hyperplanesNum, shouldReport = IO.promptUserForInput()
                    
    sub1, sub2   = hypeplanesHashing(cases, ctrls, hyperplanesNum)    #Hash SNPs with random Hyper-planes

    print "\nAfter Hashing:\n", "Cases are of size", len(sub1), "x", len(sub1[0]), "\nControls are of size", len(sub2), "x", len(sub2[0])
        
    while(1):
        decision = raw_input("Do you want to save the hashed SNPs: y:n?")
        if decision == 'y':
            IO.saveData("hashedSNPs", sub1, sub2)
            break
        if decision =='n':
            break
        
                        
    #---------------------------------------------------------------
    # Continue with KDD Algorithm
    #---------------------------------------------------------------
        
    #  Prepare Input Matrices
    original   = numpy.hstack((sub1, sub2))        # The hashed SNPs (cases, controls) concatenated.
    sub1Fliped = flipArray(sub1)                   #
    fliped     = numpy.hstack((sub1Fliped, sub2))  # As above, but with SNPs of cases flipped.
    del sub1Fliped
    
    
    print "\nStarting KDD Algorithm with Cases Flipped\n"

    initSamplingSize = len(sub1[0])/2              # The number of SNPs that will be sampled in the 1 iteration.
    results = proposePairs(cases, ctrls, original, fliped, initSamplingSize, (lowerAdaptB, upperAdaptB), shouldReport, repetitions)
    topPairStraight, totalVotesStr, avSamplingSizeStr, diffPairsEncStr = results
    
    print "Number of pairs proposed (counting duplicates)", totalVotesStr
    print "Number of Unique pairs proposed", diffPairsEncStr
    print "Average length of sampled vectors", avSamplingSizeStr
    
    #------------------------------------------------------------------------
    # KDD Algorithm, Symmetric case.
    #------------------------------------------------------------------------
            
    #prepare input matrices
    sub2Fliped = flipArray(sub2)
    fliped = numpy.hstack((sub1, sub2Fliped))  
    del (sub2Fliped, sub1, sub2)
            
    print "\nStarting KDD Algorithm with Controls Flipped\n"
    initSamplingSize = avSamplingSizeStr   #set initsamplingSize to the average of the previous runs
    results = proposePairs(cases, ctrls, original, fliped, initSamplingSize, (lowerAdaptB, upperAdaptB), shouldReport, repetitions)
    topPairSymmetric, totalVotesSym, avSamplingSizeSym, diffPairsEncSym = results
    

    print "Number of pairs proposed (counting duplicates)", totalVotesSym
    print "Number of unique pairs proposed", diffPairsEncSym
    print "Average Length of Sampled Vectors", avSamplingSizeSym
       
    print "\n\n----------------------------------------------------------------------"
    
    totalPairsProposed = totalVotesStr + totalVotesSym
    numberOfAllPairs   = len(cases) * (len(cases)-1) / float(2)

    print "We searched " + str(totalPairsProposed/numberOfAllPairs) + " of all possible SNP pairs."
    
    print "\n\n----------------------------------------------------------------------"

    toSave = raw_input("Do you want to save the proposed pairs? y:n")
    if toSave == "y":
        IO.saveData("KDD_Algo_pairs", topPairStraight, topPairSymmetric)
    
                
    ##------------------------------------------------
    ## Evaluate: comparison with Brute Force Solution
    ##------------------------------------------------
        
    solveBF = raw_input('Shall we solve it in the Brute Force way and compare the results?y:n\n')

    if solveBF == "y":        
        BFsolution = kPairsWithMaxDiffPearCorr(shouldReport, cases, ctrls)    
        straighEval = evaluateOurSolution(BFsolution, topPairStraight)
        symmetrEval = evaluateOurSolution(BFsolution, topPairSymmetric)

        top = [10, 100, 1000] 
        for tp in top:
            rec = recall_combined(straighEval, symmetrEval, tp)
            print "The KDD Algorithm found " + str(rec) + " out of the top " + str(tp) + " Best pairs."

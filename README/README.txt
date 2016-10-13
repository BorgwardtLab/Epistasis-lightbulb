Code for the paper  "Two Locus Association Mapping in Subquadratic Time", appeared in KDD 2011.  (v0.1b)
==========================================================================================================


Information
-------------------------
Author: Panagiotis Achlioptas
Affiliation: University of Crete, Heraklion / Stanford University, California
Contact: pachlioptas@gmail.com
Release date: September 18, 2012
Version: 0.1b


Installation / Usage
-------------------------
At a terminal try executing the following code:   python main.py 'casesFilePath' 'controlsFilePath' 'lowerB' 'upperB' 'repetitions' 'hyperPlanes' 'report'

Input:
'casesFilePath' and 'controlsFilePath' are the paths to the files containing the SNPs of the cases and controls respectively. 
The format of these files is expected to be the following: each line represents a single SNP. The values of the individuals corresponding to each SNP are separated by a space (" ").

For the rest of the description, let the number of SNPs be denoted as -N-.

'lowerB' and 'upperB', these are the lower and upper number of SNP pairs that you want the algorithm to propose in each repetition (hashing round). In the experiments we presented in the KDD paper, we set these values to some multiple(s) of the number of input SNPs ( e.g. 0.1 * N and 10 * N respectively ).

'repetitions' is the number of hashing rounds you want the algorithm to run. We experimentally found that setting this to the square root of the number of the input SNPs, (while maintaining the previous bounds to Theta(N)) results in sufficiently good results while reducing the overall quadratic complexity.

'hyperPlanes' is the number of hyperplanes that will be used to hash each SNP and convert it to a binary vector (its length will be equal to 'hyperPlanes'). The hashing of the input SNPs that is done with the collection of the hyperplanes, preserves the angles between the SNPs and is done only once in the beginning of the procedure.

'report' is the number of the top SNPs that you are aiming to have found at the end of the procedure. (This is the same as the size of the internal Heap that keeps the pairs with the maximum correlation difference).

For further details of how the algorithm works, please, read the Docstring of the module python.py

Concrete Example of Usage: python main.py ./cases_test ./ctrls_test 500 3000 200 200 1000


Pitfalls - Notices
-------------------------
1. Since we are searching for pairs that vary maximally between the 'Cases' and the 'Controls' according to the Pearson Correlation, each input SNP mush have non-zero variance.
2. The two input files containing the 'Cases' and the 'Controls' order the SNPs in the same manner, i.e. the same line of the two files, represent the same SNP, as measured between the different individual groups.
3. In our algorithm we use hashing in two different ways. First we use a collection of hyperplanes to convert each input SNP to a binary vector (among others: this is the way we deal with uneven sized cases/controls SNPs since after this hashing all the SNPs have the same size = number of hyperplanes used). This hashing/conversion takes place only once.
The other way involves hashing that is done repeatedly (input parameter 'repetitions' mandates how many times) and it serves in discovering SNPs that were identical in a sampled version of them.


Version history
-------------------------
Version 0.1b:
   - The initial release.


Disclaimer
-------------------------
(C) Panagiotis Achlioptas, 2010
You are free to use, modify, or redistribute this code in any way you want for non-commercial purposes. If you do so, I would appreciate it if you refer to the original author or refer to one of the papers mentioned above.


Contact
-------------------------
If you have any bugs, questions, suggestions, or modifications, please contact me:

	pachlioptas@gmail.com
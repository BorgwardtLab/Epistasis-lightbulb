# lightbulb

Efficient algorithms and GPU implementations for genome-wide epistasis screens, that is, genome-wide SNP x SNP interaction search. Please see [1] for further details.

[1] P. Achlioptas, B. Schölkopf and K. Borgwardt (2011)
**Two-locus association mapping in subquadratic runtime**,
_Proceedings of the 17th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD 2011)_ 726-734 [link](http://dl.acm.org/citation.cfm?doid=2020408.2020521)

## Usage 

The code is contained in the folder `code`. In a terminal, execute the following command:

```
python main.py 'casesFilePath' 'controlsFilePath' 'lowerB' 'upperB' 'repetitions' 'hyperPlanes' 'report'
```

where:    

* `casesFilePath` and `controlsFilePath` are the paths to the files containing the SNPs of the cases and controls respectively. The format of these files is expected to be the following: each line represents a single SNP. The values of the individuals corresponding to each SNP are separated by a space (" ").  

For the rest of the description, let the number of SNPs be denoted as `N`.  

* `lowerB` and `upperB`, these are the lower and upper number of SNP pairs that you want the algorithm to propose in each repetition (hashing round). In the experiments we presented in the KDD paper, we set these values to some multiple(s) of the number of input SNPs ( e.g. `0.1 * N` and `10 * N`, respectively).  

* `repetitions` is the number of hashing rounds you want the algorithm to run. We experimentally found that setting this to the square root of the number of the input SNPs, (while maintaining the previous bounds to `Theta(N)`) results in sufficiently good results while reducing the overall quadratic complexity.  

* `hyperPlanes` is the number of hyperplanes that will be used to hash each SNP and convert it to a binary vector (its length will be equal to `hyperPlanes`). The hashing of the input SNPs that is done with the collection of the hyperplanes, preserves the angles between the SNPs and is done only once in the beginning of the procedure.  

* `report` is the number of the top SNPs that you are aiming to have found at the end of the procedure. (This is the same as the size of the internal Heap that keeps the pairs with the maximum correlation difference).


## Sample data

There are two sample data files (simulated data, one cases and one controls) included in the `code` folder:

* `cases_test`  
* `ctrls_test`  


## Contact

For any questions or suggestions, please contact Panagiotis Achlioptas at: pachlioptas@gmail.com

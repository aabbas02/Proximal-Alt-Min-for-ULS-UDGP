# Algorithms for r-local unlabeled sensing

MATLAB code for a graph matching algorithm [1] and fast proximal alternating minimization algorithm [2] for the unlabeled sensing (ULS) problem. The algorithms are applied to the case when the unknown permutation is r-local. These codes were written by Ahmed Ali Abbasi as part of a larger research project on ULS; this project counts towards the partial fulfilment of research requirement for the PhD degree of Ahmed Ali Abbasi. The code was checked by Professor [Abiy Tasissa](http://sites.tufts.edu/atasissa/) and Professor  [Shuchin Aeron](http://www.ece.tufts.edu/~shuchin/). 

## Instructions

* The code for benchmark algorithms is contained in the folder `benchmarks`. A brief description of the benchmarks and references to the original papers are given in [2].

* To reproduce any of the results in [2], run the corerresponding file in the `figures` folder. For example, to reproduce Figure 3(a), run `main_3a.m`.

* For linux systems, replace ' \ ' in the main.m files in the `figures` folder by ' / '. These replacements need to be made in line 2 (all main files) and line 5 (main_3a.m only). 
 

## List of data files tested
* `Cropped Yale Dataset`: The dataset used for experiments in `Figure 4(a)` and `Table 1` is contained in the file `yale_compressed.mat` in the `datasets` folder. The variable `faces` of dimension 2414x96x84 contains 2414 images of dimension 96x84. 64 consecutive images are of the same face under different poses. 
* `MNIST`: The MNIST dataset used for the experiments in `Figure 4(b)` and `Table 1` is contained in the file "train set" at [URL](https://pjreddie.com/projects/mnist-in-csv/). The dataset comprises 60,000 images, each of which is represented as a 785-dimensional column vector in a matrix of dimension 60,000 x 785.  Each column vector comprises the image label (1) and the image pixeIs (2:785). In order to reproduce the results in Figure 4(b) and Table 1, download the train set and place it in the `datasets` folder.




## MATLAB files description
`make_r_local_permutation.m`: inputs are r, n; returns  an r-local permutation of size nxn with each block of size rxr. 

`lp_ls_alt_min_prox.m`: implments aternating linear program and least squares updates on variables X and P respectively. 

`lp_r_prox.m`: code for P update; the linear assignment problem is solved using the Hungarian algorithm.

`munkres.m`: input cost matrix C; return permutation P such that P = arg min <C,P>.



## References

[1]. Ahmed Ali Abbasi, Abiy Tasissa, and Shuchin Aeron. R-local unlabeled sensing: A novel graph matching approach for multiview unlabeled sensing under local permutations. IEEE Open Journal of Signal Processing, 2:309–317, 2021.
[URL](https://ieeexplore.ieee.org/document/9440727)


[2]. Ahmed Ali Abbasi, Abiy Tasissa, and Shuchin Aeron.  r-local sensing:  Improved algorithm and applications, 2021.
[URL](https://arxiv.org/abs/2110.14034)



## Feedback

Email your feedback to <a href="mailto:ahmed.abbasi@tufts.edu">Ahmed Abbasi</a>.

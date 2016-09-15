EFANNA Extreme Fast Approximate Nearest Neighbor Algorithms
============
EFANNA is a ***flexible*** and ***efficient*** library for approximate nearest neighbor search (ANN search) on large scale data. It implements the algorithms of our paper [EFANNA : Extremely Fast Approximate Nearest Neighbor Search Algorithm Based on kNN Graph](https://www.baidu.com).
EFANNA provides fast solutions on ***approximate nearest neighbor graph construction*** and ***ANN search*** based on the prebuilt graph. And it achieves best performance on million scale data.

Benchmark data set
-------
* [SIFT1M and GIST1M](http://corpus-texmex.irisa.fr/)

How To Complie
-------
Go to the root directory of EFANNA and make.

	cd efanna/
	make

How To Use
------
* An easy demo:

		cd efanna/samples/
		./sample/efanna\_sample\_kdtreeub sift\_base.fvecs  
		sift.graph sift\_query.fvecs sift.res 8 8 7 100 4   
		12 8 100 10

* Meaning of the parameters(from left to right):

	sift\_base.fvecs -- database points  
	sift.graph -- graph built by EFANNA  
	sift\_query.fvecs -- query points  
	sift.res -- ANN search results of given query  
	
	8 -- number of trees used to build the graph  
	8 -- conquer-to-depeth 
	7 -- number of iterations to build the graph 
	 
	100 -- expansion factor  
	4 -- number of iterations of ANN search  
	12 -- search-to-depth of ANN search  
	8 -- number of trees used in ANN search  
	100 -- pool size  
	10 -- required number of nearest neighbors of each query
	
See our paper or user manual for more details about the parameters.

Acknowledgment
------
Our code framework imitates [Flann](http://www.cs.ubc.ca/research/flann/) to make it scalable, and the implemnetation of NN-descent is taken from [Kgraph](http://www.kgraph.org). They proposed the NN-descent algorithm. Many thanks to them for inspiration.

What to do
-------
* Add OpenMP support for speed-up
* Add more initial algorithm choice	
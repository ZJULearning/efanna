EFANNA: an Extremely Fast Approximate Nearest Neighbor search Algorithm framework based on kNN graph
============
EFANNA is a ***flexible*** and ***efficient*** library for approximate nearest neighbor search (ANN search) on large scale data. It implements the algorithms of our paper [EFANNA : Extremely Fast Approximate Nearest Neighbor Search Algorithm Based on kNN Graph](http://arxiv.org/abs/1609.07228).    
EFANNA provides fast solutions on ***approximate nearest neighbor graph construction*** and ***ANN search*** based on the prebuilt graph. And it achieves best performance on million scale data.   
EFANNA is also flexible to adopt all kinds of hierarchical structure for initialization, such as random projection tree, hierarchical clustering tree, multi-table hashing and so on.   

What's new
-------
Interfaces change.    
AVX instructions supported.    
algorithm improved.    

Benchmark data set
-------
* [SIFT1M and GIST1M](http://corpus-texmex.irisa.fr/)

ANN search performance
------
The performance was tested without parallelism.   

![SIFT1nn](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/SIFT_1nn.png)     
![SIFT100nn](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/SIFT_100nn.png)     
![GIST1nn](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/GIST_1nn.png)     
![GIST100nn](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/GIST_100nn.png)    
Compared Algorithms:   
* [kGraph](http://www.kgraph.org)  
* [flann](http://www.cs.ubc.ca/research/flann/)   
* [IEH](http://ieeexplore.ieee.org/document/6734715/) : Fast and accurate hashing via iterative nearest neighbors expansion      
* [GNNS](https://webdocs.cs.ualberta.ca/~abbasiya/gnns.pdf) : Fast Approximate Nearest-Neighbor Search with k-Nearest Neighbor Graph     

kNN Graph Construction Performance
------
The performance was tested without parallelism.  

![SIFT1nnGraph](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/SIFT_graph.png)     
![SIFT100nnGraph](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/GIST_graph.png)    

Compared Algorithms:   
* [Kgraph](http://www.kgraph.org) (same with NN-descent)   
* [NN-expansion](https://webdocs.cs.ualberta.ca/~abbasiya/gnns.pdf) (same with GNNS)   
* [SGraph](http://ieeexplore.ieee.org/document/6247790/) : Scalable k-NN graph construction for visual descriptors  
* [FastKNN](http://link.springer.com/chapter/10.1007/978-3-642-40991-2_42) : Fast kNN Graph Construction with Locality Sensitive Hashing  
* [LargeVis](http://dl.acm.org/citation.cfm?id=2883041) : Visualizing Large-scale and High-dimensional Data    

How To Complie    
-------
Go to the root directory of EFANNA and make.    

	cd efanna/
	make

How To Use    
------
EFANNA uses a composite index to carry out ANN search, which includes an approximate kNN graph and a number of tree structures. They can be built by this library as a whole or seperately.  
  
You may build the kNN graph seperately for other use, like other graph based machine learning algorithms.  
 
 Below are some demos.   

* kNN graph building :    

		cd efanna/samples/   
		./efanna_index_buildgraph sift_base.fvecs sift.graph 8 8 8 30 25 10 10   
		  
 Meaning of the parameters(from left to right):    

	sift_base.fvecs -- database points  
	sift.graph -- graph built by EFANNA   
	
	8 -- number of trees used to build the graph (larger is more accurate but slower)   
	8 -- conquer-to-depeth(smaller is more accurate but slower)   
	8 -- number of iterations to build the graph 
	 
	30 -- L (larger is more accurate but slower, no smaller than K)  
	25 -- check (larger is more accurate but slower, no smaller than K)  
	10 -- K, for KNN graph   
	10 -- S (larger is more accurate but slower)

	
* tree building :   
    
		cd efanna/samples/
		./efanna_index_buildtrees sift_base.fvecs sift.trees 32
        
  Meaning of the parameters(from left to right):   
  
	sift_base.fvecs -- database points  
	sift.trees -- struncated KD-trees built by EFANNA  
	32 -- number of trees to build   


* index building at one time:   

		cd efanna/samples/
		./efanna_index_buildall sift_base.fvecs sift.trees sift.graph 32 8 8 200 200 100 10 8
   
  Meaning of the parameters(from left to right)   
	
	sift_base.fvecs -- database points  
	sift.trees -- struncated KD-trees built by EFANNA  
	sift.graph --
	32 -- number of trees in total for building index   
	8 -- conquer-to-depth    
	8 -- iteration number    
	200 -- L (larger is more accurate but slower, no smaller than K)   
	200 -- check (larger is more accurate but slower, no smaller than K)    
	100 -- K, for KNN graph
	10 -- S (larger is more accurate but slower)  
	8 -- 8 out of 32 trees are used for building graph

   
* ANN search

		cd efanna/samples/
		./efanna_search sift_base.fvecs sift.trees sift.graph sift_query.fvecs sift.results 16 4 1200 200 10 0
  
  Meaning of the parameters(from left to right):   
  
	sift_base.fvecs -- database points  
	sift.trees -- prebuilt struncated KD-trees used for search  
	sift.graph -- prebuilt kNN graph   
	sift_query -- sift query points  
	sift.results -- path to save ANN search results of given query   
	16 -- number of trees to use (no greater than the number of prebuilt trees)   
	4 -- number of epoches   
	1200 -- pool size factor (larger is more accurate but slower, usually 6~10 times larger than extend factor)   
	200 -- extend factor (larger is more accurate but slower)   
	10 -- required number of returned neighbors (i.e. k of k-NN)   
	0 -- searching methods (0~2, three kinds of algrothms, different performance on k-NN graph of different k, see user manual for parameter setting of other algorithms)   
 
 
* Evaluation

		cd efanna/samples/
		./evaluate sift.results sift_groundtruth.ivecs 10

  Meaning of the parameters(from left to right):   
  
	sift.results -- search results file  
	sift_groundtruth.ivecs -- ground truth file  
	10 -- evaluate the 10NN accuracy (the only first 10 points returned by the algorithm are examined, how many points are among the true 10 nearest neighbors of the query)   


See our paper or user manual for more details about the parameters and interfaces.

Output format
------
The file format of approximate kNN graph and ANN search results are the same.   
Suppose the database has N points, and numbered from 0 to N-1. You want to build an approximate kNN graph. The graph can be regarded as a N * k Matrix. The saved kNN graph binary file saves the matrix by row. The first 4 bytes of each row saves the int value of k, next 4 bytes saves the value of M and next 4 bytes saves the float value of the norm of the point. Then it follows k*4 bytes, saving the indices of the k nearest neighbors of respective point. The N rows are saved continuously without seperating characters.   

Similarly, suppose the query data has n points, numbered 0 to n-1. You want EFANNA to return k nearest neighbors for each query. The result file will save n rows like the graph file. It saves the returned indices row by row. Each row starts with 4 bytes recording value of k, and follows k*4 bytes recording neighbors' indices.  

Input of EFANNA
------
Because there is no unified format for input data, users may need to write input function to read your own data. You may imitate the input function in our sample code (sample/efanna\_efanna\_index\_buildgraph.cc) to load the data into our matrix.

To use SIMD instruction optimization, you should pay attention to the data alignment problem of SSE / AVX instruction.  

Acknowledgment
------
Our code framework imitates [Flann](http://www.cs.ubc.ca/research/flann/) to make it scalable, and the implemnetation of NN-descent is taken from [Kgraph](http://www.kgraph.org). They proposed the NN-descent algorithm. Many thanks to them for inspiration.

What to do
-------
* Add more initial algorithm choice	

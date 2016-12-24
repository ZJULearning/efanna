Search with Hash Index (and a kNN Graph)
============
We provide here the codes for search with hash index (and a kNN graph). The codes are used in our paper [A Revisit of Hashing Algorithms for Approximate Nearest Neighbor Search](http://arxiv.org/abs/1612.07545).    
Search with hash index together with a kNN graph actually is Iterative Expanding Hashing (IEH). It can be regarded as a special instance in the EFANNA framework.   
You need a hash index for search. Please use the [matlab codes](https://github.com/dengcai78/MatlabFunc/tree/master/ANNS/Hashing) of various hashing algorithms to generate the hash index.

Benchmark data set
-------
* [SIFT1M and GIST1M](http://corpus-texmex.irisa.fr/)

The performance was tested without parallelism.   

ANN search using hash index
------

![SIFT100nn](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/sift.best.time.100nn.png)     
![GIST100nn](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/gist.best.time.100nn.png)    

The number in parenthesis after the algorithm name is the optimal code length of this hashing algorithm on this dataset.

ANN search using hash index and a kNN graph (Iterative Expanding Hashing, IEH)
------

![SIFT100nn](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/sift.IEH.time.100nn.png)     
![GIST100nn](http://www.cad.zju.edu.cn/home/dengcai/Data/Hashing/gist.IEH.time.100nn.png)    

Notice that the IEH results in this figure are much better than the results reported in [EFANNA paper](https://github.com/fc731097343/efanna). It is simply because we use 16 hash tables and a 50-NN graph here, while EFANNA paper only use 1 hash table and a 10-NN graph.

How To Complie    
-------
Go to the root directory of EFANNA and make.    

	cd efanna/
	make

How To Use    
------
You need a hash index for search. Please use the [matlab codes](https://github.com/dengcai78/MatlabFunc/tree/master/ANNS/Hashing) of various hashing algorithms to generate the hash index.

* ANN search with hash index

		cd efanna/samples_hashing/
		./hashing_search sift_base.fvecs 16 LSHtableSIFT32b sift_query.fvecs LSHquerySIFT32b sift.results 16 32 8 0 10000 100

  Meaning of the parameters(from left to right):   

	sift_base.fvecs -- database points  
	LSHtableSIFT32b -- the binary code file of the database (the actual file is LSHtableSIFT32b_1, LSHtableSIFT32b_2, ...)  
	sift_query -- sift query points  
	LSHquerySIFT32b -- the binary code file of the query (the actual file is LSHquerySIFT32b_1, LSHquerySIFT32b_2, ...)  
	sift.results -- path to save ANN search results of given query   
	16 -- number of tables to use (you should provide 16 tables)   
	32 -- code length of the hash table   
	8  -- maximum radius (the algorithm will stop locate points beyond this radius number, using random points from the database to fill the initial pool)   
	0  -- code length shift (the algorithm will actually use codelength - codelengthshift bits code, in this case, 32-0=32)   
	10000 -- initial pool size factor    
	100 -- required number of returned neighbors (i.e. k of k-NN)   

To use Iterative Expanding Hashing (IEH), you need a kNN graph besides hash index. You can use [EFANNA](https://github.com/fc731097343/efanna) to build an approximate kNN graph efficiently.

* ANN search with hash index and a kNN graph (IEH)

		cd efanna/samples_hashing/
		./hashing_search sift_base.fvecs 16 LSHtableSIFT32b sift_query.fvecs LSHquerySIFT32b sift.results 16 32 8 0 200 100 sift.graph 6 100 1

  Meaning of the parameters(from left to right):   

	sift_base.fvecs -- database points  
	LSHtableSIFT32b -- the binary code file of the database (the actual file is LSHtableSIFT32b_1, LSHtableSIFT32b_2, ...)  
	sift_query -- sift query points  
	LSHquerySIFT32b -- the binary code file of the query (the actual file is LSHquerySIFT32b_1, LSHquerySIFT32b_2, ...)  
	sift.results -- path to save ANN search results of given query   
	16 -- number of tables to use (you should provide 16 tables)   
	32 -- code length of the hash table   
	8  -- maximum radius (the algorithm will stop locate points beyond this radius number, using random points from the database to fill the initial pool)   
	0  -- code length shift (the algorithm will actually use codelength - codelengthshift bits code, in this case, 32-0=32)   
  200 -- initial pool size factor    
	100 -- required number of returned neighbors (i.e. k of k-NN)   
	sift.graph -- prebuilt kNN graph   
	6 -- number of epoches   
	100 -- extend factor (larger is more accurate but slower)   
	1 -- searching methods (0~1, two kinds of algrothms, different performance on k-NN graph of different k)   


Output and Input format
------
Same as that of [EFANNA](https://github.com/fc731097343/efanna)

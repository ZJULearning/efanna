#ifndef EFANNA_DCI_INDEX_H_
#define EFANNA_DCI_INDEX_H_
#include "algorithm/base_index.hpp"
#include <fstream>
#include <time.h>
#include <string.h>
#include <random>
#include <cmath>
#include <boost/dynamic_bitset.hpp>

namespace efanna{
  struct DciIndexParams : public IndexParams
  {
	DciIndexParams(unsigned tree_total,unsigned tree_build,unsigned node_cap, int epoches = 4, int check = 25, int myL = 30, int myS = 10, int k = 10)
	{
	    init_index_type = DCI;
	    K = k;
            build_epoches = epoches;
            S = myS;
	    ValueType treev;
            treev.int_val = tree_total;
            extra_params.insert(std::make_pair("trees",treev));
            ValueType treeb;
            treeb.int_val = tree_build > 0 ? tree_build : tree_total;
            extra_params.insert(std::make_pair("treesb",treeb));
	    ValueType nodec;
            nodec.int_val = node_cap;
            extra_params.insert(std::make_pair("nodec",nodec));
	    L = myL;
            Check_K = check;
	}
  };
  template <typename DataType>
  class DciIndex : public InitIndex<DataType>
  {
  public:

    typedef InitIndex<DataType> BaseClass;
    DciIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = DciIndexParams(4,4,10)) :
      BaseClass(dataset,d,params)
    {
      std::cout<<"kdtree ub initial"<<std::endl;
      ExtraParamsMap::const_iterator it = params_.extra_params.find("trees");
      if(it != params_.extra_params.end()){
        TreeNum = (it->second).int_val;
        std::cout << "Using kdtree to build "<< TreeNum << " trees in total" << std::endl;
      }
      else{
        TreeNum = 4;
        std::cout << "Using kdtree to build "<< TreeNum << " trees in total" << std::endl;
      }
      SP.tree_num = TreeNum;

      it = params_.extra_params.find("treesb");
      if(it != params_.extra_params.end()){
        TreeNumBuild = (it->second).int_val;
        std::cout << "Building kdtree graph with "<< TreeNumBuild <<" trees"<< std::endl;
      }
      else{
        TreeNumBuild = TreeNum;
        std::cout << "Building kdtree graph with "<< TreeNumBuild <<" trees"<< std::endl;
      }

      it = params_.extra_params.find("nodec");
      if(it != params_.extra_params.end()){
        NodeC = (it->second).int_val;
        std::cout << "Building kdtree initial index with node capacity "<< NodeC  << std::endl;
      }
      else{
        NodeC = -1;
        
      }
      
    }
    void buildIndexImpl(){
      clock_t s,f;
      s = clock();
      initGraph();
      f = clock();
      std::cout << "initial graph using time: "<< (f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<< std::endl;
      //refineGraph();
    }
    void loadIndex(char* filename){}
    void saveIndex(char* filename){}
    void loadTrees(char* filename){}
    void saveTrees(char* filename){}
    void loadGraph(char* filename){
	std::ifstream in(filename,std::ios::binary);
        unsigned N;
	in.read((char*)&N, sizeof(unsigned));

        g.resize(N);
	M.resize(N);
	for(unsigned i=0; i < N; i++){
	    unsigned k,m;
	    in.read((char*)&k, sizeof(unsigned));
	    in.read((char*)&m, sizeof(unsigned));
	    M[i] = m;
	    g[i].resize(k);

	    for(unsigned j=0; j<k; j++){
		unsigned id;
		in.read((char*)&id, sizeof(unsigned));
		g[i][j].id = id;
	    }
	}
	in.close();
    }

    void saveGraph(char* filename){
	std::ofstream out(filename,std::ios::binary);
	unsigned N = g.size();
	out.write((char*)&N, sizeof(int));
	for(unsigned i=0; i < N; i++){
	    unsigned k = g[i].size();
	    unsigned m = M[i];
	    out.write((char*)&k, sizeof(unsigned));
	    out.write((char*)&m, sizeof(unsigned));
	    for(unsigned j = 0; j < k; j++){
		unsigned id = g[i][j].id;
		out.write((char*)&id, sizeof(unsigned));
	    }
	}
	out.close();
    }
    struct Node{
	unsigned StartIdx,EndIdx;
	Node* Lchild, * Rchild;

        DataType DivVal;
        ~Node() {
          if (Lchild!=NULL) Lchild->~Node();
          if (Rchild!=NULL) Rchild->~Node();
        }
    };
    void getNeighbors(size_t k, const Matrix<DataType>& query){}
    void GenDirection(std::mt19937& rng, float* d, unsigned dim){
	for (unsigned i = 0; i < dim; ++i) d[i] = rng();
	DataType n = distance_->norm(d, dim);
	n = sqrt(n);
	for (unsigned i = 0; i < dim; ++i) d[i] /= n;
    }
    void BuildTrees(Node* node, DataType* prj){
	//unsigned count = node->EndIdx - node->StartIdx;
	
    }
    void initGraph(){
	unsigned N = features_.get_rows();
	unsigned D = features_.get_cols();
	unsigned seed = 1998;
	std::mt19937 rng(seed);
        nhoods.resize(N);
        knn_graph.resize(N);
	for (auto &nhood: nhoods) {
	    //nhood.nn_new.resize(params_.S * 2);
	    nhood.pool.resize(params_.L+1);
	    nhood.radius = std::numeric_limits<float>::max();
	}
        std::vector<int> indices(N);
	LeafLists.resize(TreeNum);
	for(unsigned i = 0; i < TreeNum; i++){
	    Node* node = new Node;
            node->Lchild = NULL;
            node->Rchild = NULL;
            node->StartIdx = 0;
            node->EndIdx = N;
	    tree_roots_.push_back(node);
	    float* dir = new float[D];
	    GenDirection(rng, dir, D);
	    dirs.push_back(dir);
        }
	for(unsigned i = 0; i < TreeNum; i++){
	    std::vector<DataType> prjVals;
	    prjVals.resize(N);
	    for(unsigned j = 0; j < N; j++){
		prjVals[j] = distance_->dot(features_.get_row(j),dirs[i],D);
	    }
	    std::sort(prjVals.begin(),prjVals.end());
	    BuildTrees(tree_roots_[i], &prjVals[0]);
	}
    }
    unsigned TreeNum;
    unsigned TreeNumBuild;
    unsigned NodeC;
    std::vector<Node*> tree_roots_;
    std::vector<float*> dirs;
    std::vector<std::vector<unsigned>> LeafLists;
    USING_BASECLASS_SYMBOLS
 };


}
#endif


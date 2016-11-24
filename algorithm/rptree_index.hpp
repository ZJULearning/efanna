#ifndef EFANNA_KDTREE_UB_INDEX_H_
#define EFANNA_KDTREE_UB_INDEX_H_
#include "algorithm/base_index.hpp"
#include <fstream>
#include <time.h>
#include <string.h>
#include <random>

namespace efanna{
  struct RPtreeIndexParams : public IndexParams
  {
	RPtreeIndexParams(bool rnn_used, int tree_num_total, int merge_level, int epoches = 4, int check = 25, int myL = 30, int building_use_k = 10, int tree_num_build = 0, int myS = 10)
      {
          reverse_nn_used = rnn_used;
          init_index_type = RPTREE;
          K = building_use_k;
          build_epoches = epoches;
          S = myS;
          ValueType treev;
          treev.int_val = tree_num_total;
          extra_params.insert(std::make_pair("trees",treev));
          ValueType treeb;
          treeb.int_val = tree_num_build > 0 ? tree_num_build : tree_num_total;
          extra_params.insert(std::make_pair("treesb",treeb));
          ValueType merge_levelv;
          merge_levelv.int_val = merge_level;
          extra_params.insert(std::make_pair("ml",merge_levelv));
          L = myL;
          Check_K = check;
      }
  };
  template <typename DataType>
  class RPtreeIndex : public InitIndex<DataType>
  {
  public:
    typedef InitIndex<DataType> BaseClass;
    RPtreeIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = RPtreeIndexParams(true,4)) :
      BaseClass(dataset,d,params)
    {
	std::cout<<"rptree ub initial"<<std::endl;
      ExtraParamsMap::const_iterator it = params_.extra_params.find("trees");
      if(it != params_.extra_params.end()){
        TreeNum = (it->second).int_val;
        std::cout << "Using rptree to build "<< TreeNum << " trees in total" << std::endl;
      }
      else{
        TreeNum = 4;
        std::cout << "Using rptree to build "<< TreeNum << " trees in total" << std::endl;
      }
      SP.tree_num = TreeNum;

      it = params_.extra_params.find("treesb");
      if(it != params_.extra_params.end()){
        TreeNumBuild = (it->second).int_val;
        std::cout << "Building rptree graph with "<< TreeNumBuild <<" trees"<< std::endl;
      }
      else{
        TreeNumBuild = TreeNum;
        std::cout << "Building rptree graph with "<< TreeNumBuild <<" trees"<< std::endl;
      }

      it = params_.extra_params.find("ml");
      if(it != params_.extra_params.end()){
        ml = (it->second).int_val;
        std::cout << "Building rptree initial index with merge level "<< ml  << std::endl;
      }
      else{
        ml = -1;
        std::cout << "Building rptree initial index with max merge level "<< std::endl;
      }
      max_deepth = 0x0fffffff;
      error_flag = false;
    }
    void buildIndexImpl(){
      clock_t s,f;
      s = clock();
      initGraph();
      f = clock();
      std::cout << "initial graph using time: "<< (f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<< std::endl;

      if(error_flag){
        std::cout << "merge level deeper than tree, max merge deepth is" << max_deepth-1<<std::endl;
        return;
      }
      refineGraph();
    }
    void loadIndex(char* filename){}
    void saveIndex(char* filename){}
    void loadTrees(char* filename){}
    void saveTrees(char* filename){}
    void loadGraph(char* filename){}
    void saveGraph(char* filename){}
    struct Node
    {
        DataType* norm_vec;
	DataType offset;
        size_t StartIdx, EndIdx;
	unsigned treeid;
        Node* Lchild, * Rchild;

        ~Node() {
	  if (norm_vec!=NULL) delete[] norm_vec;
          if (Lchild!=NULL) Lchild->~Node();
          if (Rchild!=NULL) Rchild->~Node();
        }

    };

    void initGraph(){}
    int TreeNum;
    int TreeNumBuild;
    int ml;   //merge_level
    int max_deepth;
    
    //DataType* var_;
    omp_lock_t rootlock;
    bool error_flag;
    //DataType* mean_;
    std::vector<Node*> tree_roots_;
    std::vector< std::pair<Node*,size_t> > mlNodeList;
    std::vector< std::pair<Node*,size_t> > qlNodeList;
    std::vector<std::vector<unsigned>> LeafLists;
    USING_BASECLASS_SYMBOLS
  };
}
#endif

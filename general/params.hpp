#ifndef EFANNA_PARAMS_H_
#define EFANNA_PARAMS_H_
#include <map>
#include <string>

namespace efanna {
  enum init_algorithm {
      KDTREE_UB,
	  HASHING
  };

  union ValueType{
    int int_val;
    float float_val;
    char* str_pt;
  };

  typedef std::map<std::string, ValueType> ExtraParamsMap;
  struct IndexParams{

    init_algorithm init_index_type;
    size_t K;  //build knn table with nn = K
    size_t S;  //nn sets' max size
    size_t L =30;//rnn size
    size_t TNS = 10;//tree node size
    size_t Check_K = 40;
    int build_epoches;
    int extend_num; //number to extend each time
    size_t pool_size = 100;
    size_t init_num = 100;
    ExtraParamsMap extra_params;
    bool reverse_nn_used;
  };

  struct SearchParams{
    int search_init_num;
    int search_epoches;
    int search_method;
    unsigned extend_to;
    unsigned tree_num;
    unsigned search_depth;
  };

}
#endif

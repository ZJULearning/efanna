#include <efanna.hpp>
#include "handle_wrapper.hpp"
#include <map>
#include <string>
#include <sstream>

//#define char16_t LIBRARY_char16_t
#include <mex.h>
//#undef char16_t

using namespace efanna;

template<typename T>
Matrix<T> copy_matrix(const mxArray *array)
{
    int points_num = mxGetN(array);
    int dim = mxGetM(array);
    size_t mem_size = mxGetN(array)*mxGetM(array)*sizeof(T);
    void* data = malloc(mem_size);
    memcpy(data, mxGetData(array), mem_size);
    return Matrix<T>(points_num, dim, (T*)data);
}

template<typename T>
struct construct_func {
    typedef FIndex<T>* (*entry)(Matrix<T>, Distance<T>*, int, const mxArray**);
};

//TODO: default params
template<typename T>
FIndex<T>* _construct_kdtreeub(Matrix<T> dataset, Distance<T>* dist, int in_n, const mxArray *in_array[]){
    if (in_n!=8 && in_n!=9 && in_n!=2) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }

    bool rnn_used = (bool)(*mxGetPr(in_array[0]));
    int trees = (int)(*mxGetPr(in_array[1]));
    if (in_n==2) {
        mexPrintf("kdtreeub params : %d\n", trees);
        return new FIndex<T>(dataset, dist, KDTreeUbIndexParams(rnn_used, trees, 10, 10, 10, 10, 10, trees, 10));
    }

    int mlevel = (int)(*mxGetPr(in_array[2]));
    int epoches = (int)(*mxGetPr(in_array[3]));
    int L = (int)(*mxGetPr(in_array[4]));
    int check_k = (int)(*mxGetPr(in_array[5]));
    int K = (int)(*mxGetPr(in_array[6]));
    int S = (int)(*mxGetPr(in_array[7]));
    if (in_n==8) {
        mexPrintf("kdtreeub params : %d %d %d %d %d %d %d\n", trees, mlevel, epoches, L, check_k, K, S);
        return new FIndex<T>(dataset, dist, KDTreeUbIndexParams(rnn_used, trees, mlevel, epoches, check_k, L, K, trees, S));
    } else if (in_n==9) {
        int build_trees = (int)(*mxGetPr(in_array[8]));
        mexPrintf("kdtreeub params : %d %d %d %d %d %d %d %d\n", trees, mlevel, epoches, L, check_k, K, S, build_trees);
        return new FIndex<T>(dataset, dist, KDTreeUbIndexParams(rnn_used, trees, mlevel, epoches, check_k, L, K, build_trees, S));
    }   
    return NULL; // ERROR if this line is reached
}

/*
template<typename T>
FIndex<T>* _construct_nndescent(Matrix<T> dataset, Distance<T>* dist, int in_n, const mxArray *in_array[]){
    if (in_n!=4) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }
    bool rnn_used = (bool)(*mxGetPr(in_array[0]));
    int epoches = (int)(*mxGetPr(in_array[1]));
    int K = (int)(*mxGetPr(in_array[2]));
    int L = (int)(*mxGetPr(in_array[3]));
    mexPrintf("nndescent params : %d %d %d %d\n", rnn_used, epoches, K, L);
    return new FIndex<T>(dataset, dist, RandomIndexParams(rnn_used, epoches, K, L));
}

template<typename T>
FIndex<T>* _construct_nnexp(Matrix<T> dataset, Distance<T>* dist, int in_n, const mxArray *in_array[]){
    if (in_n!=2) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }
    int epoches = (int)(*mxGetPr(in_array[0]));
    int extend = (int)(*mxGetPr(in_array[1]));
    mexPrintf("nnexp params : %d %d\n", epoches, extend);
    return new FIndex<T>(dataset, dist, NNexpIndexParams(epoches, extend));
}
*/

template<typename T>
void _construct(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n<3) || (!mxIsChar(in_array[1])) || (!mxIsChar(in_array[2]))) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }
    if (out_n!=1) {
        mexErrMsgTxt("Incorrect number of output arguments");
    }

    std::map<std::string, Distance<T>* > dist_table = {
        {"l2", new L2DistanceAVX<T>() }
    };

    std::map<std::string, typename construct_func<T>::entry> index_table = {
        {"kdtreeub", &_construct_kdtreeub<T> },
        {"nndescent", &_construct_kdtreeub<T> },
        {"nnexp", &_construct_kdtreeub<T> }
    };
/*
    std::map<std::string, typename construct_func<T>::entry> index_table = {
        {"kdtreeub", &_construct_kdtreeub<T> },
        {"nndescent", &_construct_nndescent<T> },
        {"nnexp", &_construct_nnexp<T> }
    };
*/

    std::string index_name = mxArrayToString(in_array[1]);
    if (index_table.find(index_name)==index_table.end()) {
        mexErrMsgTxt("Error: bad distance selector: wrong distance name.");
    }
    std::string dist_name = mxArrayToString(in_array[2]);
    if (dist_table.find(dist_name)==dist_table.end()) {
        mexErrMsgTxt("Error: bad distance selector: wrong distance name.");
    }

    //TODO: type check
    Matrix<T> dataset = copy_matrix<T>(in_array[0]);

    FIndex<T>* result = (*index_table[index_name])(dataset, dist_table[dist_name], in_n-3, in_array+3);
    out_array[0] = handle2mat<FIndex<T> >(result);
    //mxFree(index_name);
    //mxFree(dist_name);
}

template<typename T>
void _destruct(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if (in_n !=1) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }
    delete_wrapper<FIndex<T> >(in_array[0]);
}

template<typename T>
void _get_graph_mat(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    //TODO: verify there're no bugs
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    size_t nrows = handle->getGraphSize();
    std::map<unsigned, std::vector<unsigned>*> col_graph;
    int maxnnz = 0;
    for (unsigned i = 0; i < nrows; i++) {
        std::vector<unsigned> nodes = handle->getGraphRow(i);
        maxnnz += nodes.size();
        for (unsigned x : nodes) {
            if (col_graph.find(x) == col_graph.end()) {
                col_graph[x] = new std::vector<unsigned>();
            }
            col_graph[x]->push_back(i);
        }
    }
    
    //the type must be double*
    double* pr = (double *)mxCalloc(maxnnz, sizeof(double));
    mwIndex* ir = (size_t *)mxCalloc(maxnnz, sizeof(mwIndex));
    mwIndex* jc = (size_t *)mxCalloc(nrows + 1, sizeof(mwIndex));
    int nfilled = -1;
    int njc = 0;
    jc[0] = 0;
    for (unsigned i = 0; i < nrows; i++) {
        if (col_graph.find(i) == col_graph.end()) {
            njc ++;
            jc[njc] = jc[njc-1];
            continue;
        }
        for (unsigned x: *col_graph[i]) {
            nfilled ++;
            pr[nfilled] = 1;
            ir[nfilled] = x;
        }
        njc ++;
        jc[njc] = jc[njc-1] + col_graph[i]->size();
        delete col_graph[i];
    }

    out_array[0] = mxCreateSparse(nrows, nrows, maxnnz, mxREAL);
    mxSetPr(out_array[0], pr);
    mxSetIr(out_array[0], ir);
    mxSetJc(out_array[0], jc);
}

template<typename T>
void _build_index(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if (in_n != 1) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }

    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->buildIndex();
    _get_graph_mat<T>(out_n, out_array, in_n, in_array);
}

template<typename T>
void _build_trees(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if (in_n != 1) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }

    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->buildTrees();
}

template<typename T>
void _load_index(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n != 2) || (!mxIsChar(in_array[1]))) {
        mexErrMsgTxt("Error: bad path name");
    }
    char* path = mxArrayToString(in_array[1]);
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->loadIndex(path);
    mxFree(path);
}

template<typename T>
void _save_index(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n != 2) || (!mxIsChar(in_array[1]))) {
        mexErrMsgTxt("Error: bad path name");
    }
    char* path = mxArrayToString(in_array[1]);
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->saveIndex(path);
    mxFree(path);
}

template<typename T>
void _load_graph(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n != 2) || (!mxIsChar(in_array[1]))) {
        mexErrMsgTxt("Error: bad path name");
    }
    char* path = mxArrayToString(in_array[1]);
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->loadGraph(path);
    _get_graph_mat<T>(out_n, out_array, in_n, in_array);
    mxFree(path);
}

template<typename T>
void _save_graph(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n != 2) || (!mxIsChar(in_array[1]))) {
        mexErrMsgTxt("Error: bad path name");
    }
    char* path = mxArrayToString(in_array[1]);
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->saveGraph(path);
    mxFree(path);
}

template<typename T>
void _load_trees(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n != 2) || (!mxIsChar(in_array[1]))) {
        mexErrMsgTxt("Error: bad path name");
    }
    char* path = mxArrayToString(in_array[1]);
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->loadTrees(path);
    mxFree(path);
}

template<typename T>
void _save_trees(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n != 2) || (!mxIsChar(in_array[1]))) {
        mexErrMsgTxt("Error: bad path name");
    }
    char* path = mxArrayToString(in_array[1]);
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->saveTrees(path);
    mxFree(path);
}

template<typename T>
void _set_search_params_kdtreeub(FIndex<T>* handle, int in_n, const mxArray *in_array[]) {
    int search_trees = (int)(*mxGetPr(in_array[0]));
    int search_epoc = (int)(*mxGetPr(in_array[1]));
    int search_extend = (int)(*mxGetPr(in_array[2]));
    int poolsz = (int)(*mxGetPr(in_array[3]));
    int search_method = (int)(*mxGetPr(in_array[4]));
    handle->setSearchParams(search_epoc, poolsz, search_extend, search_trees, search_method);
}

template<typename T>
void _set_search_params_nndescent(FIndex<T>* handle, int in_n, const mxArray *in_array[]) {
    mexPrintf("No param needs to be set.\n");
}

template<typename T>
void _set_search_params_nnexp(FIndex<T>* handle, int in_n, const mxArray *in_array[]) {
    mexPrintf("No param needs to be set.\n");
}

template<typename T>
void _set_search_params(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n < 2) || (!mxIsChar(in_array[1]))) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }

    //XXX: why this does not need template type for the function?
    typedef void (*funcp)(FIndex<T>*, int, const mxArray**);
    std::map<std::string, funcp> index_table = {
        {"kdtreeub", &_set_search_params_kdtreeub},
        {"nndescent", &_set_search_params_nndescent},
        {"nnexp", &_set_search_params_nnexp}
    };
    
    std::string index_name = mxArrayToString(in_array[1]);
    if (index_table.find(index_name)==index_table.end()) {
        mexErrMsgTxt("Error: bad distance selector: wrong distance name.");
    }
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    (*index_table[index_name])(handle, in_n-2, in_array+2);
    //mxFree(index_name);
}

template<typename T>
void _knn_search(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    int k = (int)(*mxGetPr(in_array[1]));
    //TODO: type check
    Matrix<T> query = copy_matrix<T>(in_array[2]);
    handle->knnSearch(k, query);
}

template<typename T>
void _save_result(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    if ((in_n != 2) || (!mxIsChar(in_array[1]))) {
        mexErrMsgTxt("Error: bad path name");
    }
    char* path = mxArrayToString(in_array[1]);
    FIndex<T>* handle = mat2handle<FIndex<T> >(in_array[0]);
    handle->saveResults(path);
    mxFree(path);
}

template<typename T>
struct mexfunc {
    typedef void (*entry)(int, mxArray**, int, const mxArray**);
};

void mexFunction(int out_n, mxArray* out_array[], int in_n, const mxArray *in_array[]) {
    
    std::map<std::string, mexfunc<float>::entry> ftable = {
        {"new", &_construct<float>},
        {"destruct", &_destruct<float>},
        {"build_index", &_build_index<float>},
        {"build_trees", &_build_trees<float>},
        {"load_index", &_load_index<float>},
        {"save_index", &_save_index<float>},
        {"load_graph", &_load_graph<float>},
        {"save_graph", &_save_graph<float>},
        {"load_trees", &_load_trees<float>},
        {"save_trees", &_save_trees<float>},
        {"set_search_params", &_set_search_params<float>}, 
        {"knn_search", &_knn_search<float>}, 
        {"save_result", &_save_result<float>} 
//      {"set_distance_type", &_set_distance_type} //TODO: only L2 distance is provided in general/distance.hpp now
    };

    if ((in_n<=1) || (!mxIsChar(in_array[0]))) {
        mexErrMsgTxt("Error: bad function selector: wrong type of function name.");
    }
    std::string func_name = mxArrayToString(in_array[0]);
    if (ftable.find(func_name)==ftable.end()) {
        mexErrMsgTxt("Error: bad function selector: wrong function name.");
    }
    mexPrintf("Running Function : %s ...\n", func_name.c_str());
    (*ftable[func_name])(out_n, out_array, in_n-1, in_array+1);
    //mxFree(func_name);
}

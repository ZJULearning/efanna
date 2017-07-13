#include <efanna.hpp>
#include <iostream>
#include <fstream>
#include <ctime>
#include <malloc.h>
#include <boost/timer/timer.hpp>

using namespace boost;
using namespace efanna;
using namespace std;
void load_data(char* filename, float*& data, size_t& num,int& dim){// load data with sift10K pattern
  ifstream in(filename, ios::binary);
  if(!in.is_open()){cout<<"open file error"<<endl;exit(-1);}
  in.read((char*)&dim,4);
  cout<<"data dimension: "<<dim<<endl;
  in.seekg(0,ios::end);
  ios::pos_type ss = in.tellg();
  size_t fsize = (size_t)ss;
  num = fsize / (dim+1) / 4;
  int cols = (dim + 7)/8*8;
  data = (float*)memalign(KGRAPH_MATRIX_ALIGN, num * cols * sizeof(float));
if(dim!=cols)cout<<"data align to dimension "<<cols<<" for avx2 inst"<<endl;

  in.seekg(0,ios::beg);
  for(size_t i = 0; i < num; i++){
    in.seekg(4,ios::cur);
    in.read((char*)(data+i*cols),dim*4);
  }
  in.close();
}
int main(int argc, char** argv){
  if(argc!=11 && argc!=12){cout<< argv[0] << " data_file tree_index graph_index query_file result_file tree epoch initsz extend querNN search_method(optional)" <<endl; exit(-1);}

  float* data_load = NULL;
  float* query_load = NULL;
  size_t points_num;
  int dim;
  load_data(argv[1], data_load, points_num,dim);
  size_t q_num;
  int qdim;
  load_data(argv[4], query_load, q_num,qdim);
  Matrix<float> dataset(points_num,dim,data_load);
  Matrix<float> query(q_num,qdim,query_load);

  FIndex<float> index(dataset, new L2DistanceAVX<float>(), efanna::KDTreeUbIndexParams(true, 8 ,8 ,10,25,30,10));
  index.loadTrees(argv[2]);
  index.loadGraph(argv[3]);

  int search_trees = atoi(argv[6]);
  int search_epoc = atoi(argv[7]);
  int poolsz = atoi(argv[8]);
  int search_extend = atoi(argv[9]);
  int search_method = argc == 12 ? atoi(argv[11]) : 0;
  index.setSearchParams(search_epoc, poolsz, search_extend, search_trees, search_method);


//clock_t s,f;
boost::timer::auto_cpu_timer timer;
//s=clock();
  index.knnSearch(atoi(argv[10])/*query nn*/,query);
//f=clock();

//cout<<(f-s)*1.0/CLOCKS_PER_SEC/8<<endl;
cout<<timer.elapsed().wall / 1e9<<endl;
  index.saveResults(argv[5]);
  return 0;
}

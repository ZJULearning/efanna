#include <efanna.hpp>
#include <iostream>
#include <fstream>
#include <ctime>
#include <malloc.h>

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
  int cols = (dim + 3)/4*4;
  data = (float*)memalign(KGRAPH_MATRIX_ALIGN, num * cols * sizeof(float));
if(dim!=cols)cout<<"data align to dimension "<<cols<<" for sse2 inst"<<endl;

  in.seekg(0,ios::beg);
  for(size_t i = 0; i < num; i++){
    in.seekg(4,ios::cur);
    in.read((char*)(data+i*cols),dim*4);
  }
  in.close();
}
int main(int argc, char** argv){
  if(argc!=14){cout<<"arguments error"<<endl; exit(-1);}
  float* data_load = NULL;
  float* query_load = NULL;
  size_t points_num;
  int dim;
  load_data(argv[1], data_load, points_num,dim);
  size_t q_num;
  int qdim;
  load_data(argv[3], query_load, q_num,qdim);
  Matrix<float> dataset(points_num,dim,data_load);
  Matrix<float> query(q_num,qdim,query_load);

  unsigned int trees = atoi(argv[5]);
  int mlevel = atoi(argv[6]);
  unsigned int epochs = atoi(argv[7]);
  FIndex<float> index(dataset, new L2Distance<float>(), efanna::KDTreeUbIndexParams(true, trees ,mlevel ,epochs,25,30,10));
  clock_t s,f;
  s = clock();
  index.buildIndex();

  f = clock();
  cout<<"Index building time : "<<(f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<<endl;
  //index.saveIndex(argv[2]);
  //index.loadIndex(argv[2]);
  //index.saveTrees(argv[2]);
  //index.loadTrees(argv[2]);
  //index.loadGraph("../sift/sift_bf.index");
  index.saveGraph(argv[2]);
  int search_epoc = atoi(argv[9]);
  int search_extend = atoi(argv[8]);
  int search_trees = atoi(argv[11]);
  int search_lv = atoi(argv[10]);
  int poolsz = atoi(argv[12]);
  index.setSearchParams(search_epoc, poolsz, search_extend, search_trees,search_lv);
  s = clock();
  index.knnSearch(atoi(argv[13])/*query nn*/,query);
  f = clock();
  cout<<"Query searching time : "<<(f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<<endl;
  index.saveResults(argv[4]);
  return 0;
}

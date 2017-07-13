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
  if(argc!=4){cout<< argv[0] << " data_file save_trees_file trees" <<endl; exit(-1);}

  float* data_load = NULL;
  //float* query_load = NULL;
  size_t points_num;
  int dim;
  load_data(argv[1], data_load, points_num,dim);
  Matrix<float> dataset(points_num,dim,data_load);

  unsigned int trees = atoi(argv[3]);

  FIndex<float> index(dataset, new L2DistanceAVX<float>(), efanna::KDTreeUbIndexParams(true, trees ,8 ,8,25,30, 10, 8, 10));
  clock_t s,f;
  s = clock();
  index.buildTrees();

  f = clock();
  cout<<"Index building time : "<<(f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<<endl;
  index.saveTrees(argv[2]);
  return 0;
}

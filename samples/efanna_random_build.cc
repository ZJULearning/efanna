#include <efanna.hpp>
#include <iostream>
#include <fstream>
#include <ctime>

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
  data = new float[num*dim];

  in.seekg(0,ios::beg);
  for(size_t i = 0; i < num; i++){
    in.seekg(4,ios::cur);
    in.read((char*)(data+i*dim),dim*4);
  }
  in.close();
}
int main(int argc, char** argv){
  if(argc!=8){cout<<"arguments error"<<endl; exit(-1);}
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
  int epoch = atoi(argv[5]);
  int K = atoi(argv[6]);
  int L = atoi(argv[7]);
  FIndex<float> index(dataset, new L2Distance<float>(), efanna::RandomIndexParams(true, epoch, K, L));
  clock_t s,f;
  s = clock();
  index.buildIndex();
  f = clock();
  cout<<"Index building time : "<<(f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<<endl;
  index.saveIndex(argv[2]);
  index.setSearchParams(5, 100, 25);
  s = clock();
  index.knnSearch(100,query);
  f = clock();
  cout<<"Query searching time : "<<(f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<<endl;
  index.saveResults(argv[4]);
  return 0;
}

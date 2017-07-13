#include <efanna.hpp>
#include <iostream>
#include <fstream>
#include <ctime>
#include <malloc.h>
#include <boost/timer/timer.hpp>

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
  if(argc!=10){cout<< argv[0] << " data_file save_graph_file trees level epoch L K kNN S" <<endl; exit(-1);}

  float* data_load = NULL;
  //float* query_load = NULL;
  size_t points_num;
  int dim;
  load_data(argv[1], data_load, points_num,dim);
  //size_t q_num;
  //int qdim;
  //load_data(argv[3], query_load, q_num,qdim);
  Matrix<float> dataset(points_num,dim,data_load);
  //Matrix<float> query(q_num,qdim,query_load);

  unsigned int trees = atoi(argv[3]);
  int mlevel = atoi(argv[4]);
  unsigned int epochs = atoi(argv[5]);
  int L = atoi(argv[6]);
  int checkK = atoi(argv[7]);
  int kNN = atoi(argv[8]);
  int S = atoi(argv[9]);

  //srand(time(NULL));
  FIndex<float> index(dataset, new L2DistanceAVX<float>(), efanna::KDTreeUbIndexParams(true, trees ,mlevel ,epochs,checkK,L, kNN, trees, S));
boost::timer::auto_cpu_timer timer;

  index.buildIndex();
cout<<timer.elapsed().wall / 1e9<<endl;
  index.saveGraph(argv[2]);
  return 0;
}

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
	if(argc!=13){cout<< argv[0] << " data_file UpperBits BaseCodeFile graph_index query_file QueryCodeFile result_file radius epoch initsz extend querNN" <<endl; exit(-1);}

  float* data_load = NULL;
  float* query_load = NULL;
  size_t points_num;
  int dim;
  load_data(argv[1], data_load, points_num,dim);
  size_t q_num;
  int qdim;
  load_data(argv[5], query_load, q_num,qdim);
  Matrix<float> dataset(points_num,dim,data_load);
  Matrix<float> query(q_num,qdim,query_load);

  int UpperBits = atoi(argv[2]);
  int radius = atoi(argv[8]);
  FIndex<float> index(dataset, new L2Distance<float>(), efanna::IEHIndexParams(UpperBits,radius,argv[3],argv[6]));
  index.loadGraph(argv[4]);

  clock_t s,f;
  s = clock();
  index.buildIndex();
  f = clock();
  cout<<"Index building time : "<<(f-s)/CLOCKS_PER_SEC<<" seconds"<<endl;

  int search_epoc = atoi(argv[9]);
  int poolsz = atoi(argv[10]);
  int search_extend = atoi(argv[11]);
  index.setSearchParams(search_epoc, poolsz, search_extend);

  s = clock();
  index.knnSearch(atoi(argv[12])/*query nn*/,query);
  f = clock();
  cout<<"Query searching time : "<<(f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<<endl;
  index.saveResults(argv[7]);
  return 0;
}

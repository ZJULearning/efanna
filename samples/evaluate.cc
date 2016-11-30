#include <iostream>
#include <fstream>
#include <set>

using namespace std;
void load_data(char* filename, int*& data, int& num,int& dim){// load data with sift10K pattern
  ifstream in(filename, ios::binary);
  if(!in.is_open()){cout<<"open file error"<<endl;return;}
  in.read((char*)&dim,4);
  in.seekg(0,ios::end);
  ios::pos_type ss = in.tellg();
  int fsize = (int)ss;
  num = fsize / (dim+1) / 4;
  data = new int[num*dim];

  in.seekg(0,ios::beg);
  for(int i = 0; i < num; i++){
    in.seekg(4,ios::cur);
    in.read((char*)(data+i*dim),dim*4);
  }
  in.close();
}

int main(int argc, char** argv){
  if(argc!=4){cout<< argv[0] << " resultfile gound_truth kNN" <<endl; exit(-1);}
  int *gt = NULL;
  int *gt_ = NULL;
  int p1,p2;
  int dim1,dim2;
  load_data(argv[1], gt,  p1,dim1);
  load_data(argv[2], gt_, p2,dim2);

  if(p1 != p2){cout<< "result file and groundtruth don't match" <<endl; exit(-1);}

  int kNN = atoi(argv[3]);

  if(kNN > dim1){cout<< "result file not enough for k="<<kNN <<endl; exit(-1);}

  set<int> result;
  for(int i=0; i < p1; i++){
    result.clear();
    for(int j=0; j<kNN;j++){
      result.insert(gt[i*dim1+j]);
    }
    if (result.size()< (unsigned int) kNN)
      cout<< "query " << i <<" result index not unique!"<< endl;

  }

  int cnt=0;
  for(int i = 0; i < p1; i++){
    for(int j=0; j< kNN; j++){
      int k=0;
      for(; k<kNN; k++){
        if(gt[i*dim1+j]==gt_[i*dim2+k])break;
      }
      if(k==kNN)cnt++;
    }
  }
  cout<<kNN <<"NN accuracy: "<<1-(float)cnt/(p1*kNN)<<endl;


  return 0;

}

#ifndef EFANNA_HASHING_INDEX_H_
#define EFANNA_HASHING_INDEX_H_

#include "algorithm/base_index.hpp"
#include <boost/dynamic_bitset.hpp>
#include <time.h>
//for Debug
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <map>
#include <sstream>
#include <set>
//#define MAX_RADIUS 6


namespace efanna{
struct HASHINGIndexParams : public IndexParams
{
    HASHINGIndexParams(int codelen, int TableNum,int UpperBits, int HashRadius, char*& BaseCodeFile, char*& QueryCodeFile, int codelenShift = 0)
    {
      init_index_type = HASHING;
      ValueType len;
      len.int_val = codelen;
      extra_params.insert(std::make_pair("codelen",len));
      ValueType nTab;
      nTab.int_val = TableNum;
      extra_params.insert(std::make_pair("tablenum",nTab));
      ValueType upb;
      upb.int_val = UpperBits;
      extra_params.insert(std::make_pair("upbits",upb));
      ValueType radius;
      radius.int_val = HashRadius;
      extra_params.insert(std::make_pair("radius",radius));
      ValueType bcf;
      bcf.str_pt = BaseCodeFile;
      extra_params.insert(std::make_pair("bcfile",bcf));
      ValueType qcf;
      qcf.str_pt = QueryCodeFile;
      extra_params.insert(std::make_pair("qcfile",qcf));
      ValueType lenShift;
      lenShift.int_val = codelenShift;
      extra_params.insert(std::make_pair("lenshift",lenShift));
    }
};

template <typename DataType>
class HASHINGIndex : public InitIndex<DataType>
{
  public:

    typedef InitIndex<DataType> BaseClass;

    typedef std::vector<unsigned int> Codes;
    typedef std::unordered_map<unsigned int, std::vector<unsigned int> > HashBucket;
    typedef std::vector<HashBucket> HashTable;

    typedef std::vector<unsigned long> Codes64;
    typedef std::unordered_map<unsigned long, std::vector<unsigned int> > HashBucket64;
    typedef std::vector<HashBucket64> HashTable64;


    HASHINGIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = HASHINGIndexParams(0,NULL,NULL)) :
      BaseClass(dataset,d,params)
    {
      std::cout<<"HASHING initial, max code length : 64" <<std::endl;

      ExtraParamsMap::const_iterator it = params_.extra_params.find("codelen");
      if(it != params_.extra_params.end()){
        codelength = (it->second).int_val;
        std::cout << "use  "<<codelength<< " bit code"<< std::endl;
      }
      else{
        std::cout << "error: no code length setting" << std::endl;
      }

      it = params_.extra_params.find("lenshift");
      if(it != params_.extra_params.end()){
        codelengthshift = (it->second).int_val;
        int actuallen = codelength - codelengthshift;
        if(actuallen > 0){
          std::cout << "Actually use  "<< actuallen<< " bit code"<< std::endl;
        }else{
          std::cout << "lenShift error: could not be larger than the code length!  "<<  std::endl;
        }
      }
      else{
        codelengthshift = 0;
      }

      it = params_.extra_params.find("tablenum");
      if(it != params_.extra_params.end()){
        tablenum = (it->second).int_val;
        std::cout << "use  "<<tablenum<< " hashtables"<< std::endl;
      }
      else{
        std::cout << "error: no table number setting" << std::endl;
      }

      it = params_.extra_params.find("upbits");
      if(it != params_.extra_params.end()){
        upbits = (it->second).int_val;
        std::cout << "use upper "<<upbits<< " bits as first level index of hashtable"<< std::endl;
        std::cout << "use lower "<<codelength - codelengthshift - upbits<< " bits as second level index of hashtable"<< std::endl;
      }
      else{
        std::cout << "error: no upper bits number setting" << std::endl;
      }

      if(upbits >= codelength-codelengthshift){
        std::cout << "upbits should be smaller than the actual codelength!" << std::endl;
        return;
      }


      int actuallen = codelength - codelengthshift;
      it = params_.extra_params.find("radius");
      if(it != params_.extra_params.end()){
        radius = (it->second).int_val;
        if(actuallen<=32){
          if(radius > 13){
            std::cout << "radius greater than 13 not supported yet!" << std::endl;
            radius = 13;
          }
        }else if(actuallen<=36){
          if(radius > 11){
            std::cout << "radius greater than 11 not supported yet!" << std::endl;
            radius = 11;
          }
        }else if(actuallen<=40){
          if(radius > 10){
            std::cout << "radius greater than 10 not supported yet!" << std::endl;
            radius = 10;
          }
        }else if(actuallen<=48){
          if(radius > 9){
            std::cout << "radius greater than 9 not supported yet!" << std::endl;
            radius = 9;
          }
        }else if(actuallen<=60){
          if(radius > 8){
            std::cout << "radius greater than 8 not supported yet!" << std::endl;
            radius = 8;
          }
        }else{ //actuallen<=64
          if(radius > 7){
            std::cout << "radius greater than 7 not supported yet!" << std::endl;
            radius = 7;
          }
        }
        std::cout << "search hamming radius "<<radius<< std::endl;
      }
      else{
        std::cout << "error: no radius number setting" << std::endl;
      }

      it = params_.extra_params.find("bcfile");
      if(it != params_.extra_params.end()){
        char* fpath = (it->second).str_pt;
        std::string str(fpath);
        std::cout << "Loading base code from " << str << std::endl;

        if (codelength <= 32 ){
          LoadCode32(fpath, BaseCode);
        }else if(codelength <= 64 ){
          LoadCode64(fpath, BaseCode64);
        }else{
          std::cout<<"code length not supported yet!"<<std::endl;
        }

        std::cout << "code length is "<<codelength<<std::endl;
      }
      else{
        std::cout << "error: no base code file" << std::endl;
      }

      it = params_.extra_params.find("qcfile");
      if(it != params_.extra_params.end()){
        char* fpath = (it->second).str_pt;
        std::string str(fpath);
        std::cout << "Loading query code from " << str << std::endl;

        if (codelength <= 32 ){
          LoadCode32(fpath, QueryCode);
        }else if(codelength <= 64 ){
          LoadCode64(fpath, QueryCode64);
        }else{
          std::cout<<"code length not supported yet!"<<std::endl;
        }

        std::cout << "code length is "<<codelength<<std::endl;
      }
      else{
        std::cout << "error: no query code file" << std::endl;
      }

    }



    /**
     * Builds the index
     */

    void LoadCode32(char* filename, std::vector<Codes>& baseAll){
      if (tablenum < 1){
        std::cout<<"Total hash table num error! "<<std::endl;
      }

      int actuallen = codelength-codelengthshift;

      unsigned int  maxValue = 1;
      for(int i=1;i<actuallen;i++){
        maxValue = maxValue << 1;
        maxValue ++;
      }


      std::stringstream ss;
      for(int j = 0; j < tablenum; j++){

        ss << filename << "_" << j+1 ;
        std::string codeFile;
        ss >> codeFile;
        ss.clear();

        std::ifstream in(codeFile.c_str(), std::ios::binary);
        if(!in.is_open()){std::cout<<"open file " << filename <<" error"<< std::endl;return;}

        int codeNum;
        in.read((char*)&codeNum,4);
        if (codeNum != 1){
          std::cout<<"Codefile  "<< j << " error!"<<std::endl;
        }

        in.read((char*)&codelength,4);
        //std::cout<<"codelength: "<<codelength<<std::endl;

        int num;
        in.read((char*)&num,4);
        //std::cout<<"ponit num: "<<num<<std::endl;

        Codes base;
        for(int i = 0; i < num; i++){
          unsigned int codetmp;
          in.read((char*)&codetmp,4);
          codetmp = codetmp >> codelengthshift;
          if (codetmp > maxValue){
            std::cout<<"codetmp: "<< codetmp <<std::endl;
            std::cout<<"codelengthshift: "<<codelengthshift<<std::endl;
            std::cout<<"codefile  "<< codeFile << " error! Exceed maximum value"<<std::endl;
            in.close();
            return;
          }
          base.push_back(codetmp);
        }
        baseAll.push_back(base);

        in.close();
      }
    }

    void LoadCode64(char* filename, std::vector<Codes64>& baseAll){
      if (tablenum < 1){
        std::cout<<"Total hash table num error! "<<std::endl;
      }

      int actuallen = codelength-codelengthshift;

      unsigned long  maxValue = 1;
      for(int i=1;i<actuallen;i++){
        maxValue = maxValue << 1;
        maxValue ++;
      }


      std::stringstream ss;
      for(int j = 0; j < tablenum; j++){

        ss << filename << "_" << j+1 ;
        std::string codeFile;
        ss >> codeFile;
        ss.clear();

        std::ifstream in(codeFile.c_str(), std::ios::binary);
        if(!in.is_open()){std::cout<<"open file " << filename <<" error"<< std::endl;return;}

        int codeNum;
        in.read((char*)&codeNum,4);
        if (codeNum != 1){
          std::cout<<"Codefile  "<< j << " error!"<<std::endl;
        }

        in.read((char*)&codelength,4);
        //std::cout<<"codelength: "<<codelength<<std::endl;

        int num;
        in.read((char*)&num,4);
        //std::cout<<"ponit num: "<<num<<std::endl;

        Codes64 base;
        for(int i = 0; i < num; i++){
          unsigned long codetmp;
          in.read((char*)&codetmp,8);
          codetmp = codetmp >> codelengthshift;
          if (codetmp > maxValue){
            std::cout<<"codetmp: "<< codetmp <<std::endl;
            std::cout<<"codelengthshift: "<<codelengthshift<<std::endl;
            std::cout<<"codefile  "<< codeFile << " error! Exceed maximum value"<<std::endl;
            in.close();
            return;
          }

          base.push_back(codetmp);
        }
        baseAll.push_back(base);

        in.close();
      }
    }

    void BuildHashTable32(int upbits, int lowbits, std::vector<Codes>& baseAll ,std::vector<HashTable>& tbAll){

      for(size_t h=0; h < baseAll.size(); h++){
        Codes& base = baseAll[h];

        HashTable tb;
        for(int i = 0; i < (1 << upbits); i++){
          HashBucket emptyBucket;
          tb.push_back(emptyBucket);
        }

        for(size_t i = 0; i < base.size(); i ++){
          unsigned int idx1 = base[i] >> lowbits;
          unsigned int idx2 = base[i] - (idx1 << lowbits);
          if(tb[idx1].find(idx2) != tb[idx1].end()){
            tb[idx1][idx2].push_back(i);
          }else{
            std::vector<unsigned int> v;
            v.push_back(i);
            tb[idx1].insert(make_pair(idx2,v));
          }
        }
        tbAll.push_back(tb);
      }
    }

    void generateMask32(){
      //i = 0 means the origin code
      HammingBallMask.push_back(0);
      HammingRadius.push_back(HammingBallMask.size());

      if(radius>0){
        //radius 1
        for(int i = 0; i < codelength; i++){
          unsigned int mask = 1 << i;
          HammingBallMask.push_back(mask);
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>1){
        //radius 2
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            unsigned int mask = (1<<i) | (1<<j);
            HammingBallMask.push_back(mask);
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>2){
        //radius 3
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              unsigned int mask = (1<<i) | (1<<j) | (1<<k);
              HammingBallMask.push_back(mask);
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>3){
        //radius 4
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a);
                HammingBallMask.push_back(mask);
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>4){
        //radius 5
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b);
                  HammingBallMask.push_back(mask);
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>5){
        //radius 6
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b)| (1<<c);
                    HammingBallMask.push_back(mask);
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>6){
        //radius 7
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b)| (1<<c)| (1<<d);
                      HammingBallMask.push_back(mask);
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>7){
        //radius 8
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b)| (1<<c)| (1<<d)| (1<<e);
                        HammingBallMask.push_back(mask);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>8){
        //radius 9
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        for(int f = e+1; f < codelength; f++){
                          unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b)| (1<<c)| (1<<d)| (1<<e)| (1<<f);
                          HammingBallMask.push_back(mask);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>9){
        //radius 10
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        for(int f = e+1; f < codelength; f++){
                          for(int g = f+1; g < codelength; g++){
                            unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b)| (1<<c)| (1<<d)| (1<<e)| (1<<f)| (1<<g);
                            HammingBallMask.push_back(mask);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>10){
        //radius 11
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        for(int f = e+1; f < codelength; f++){
                          for(int g = f+1; g < codelength; g++){
                            for(int h = g+1; h < codelength; h++){
                              unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b)| (1<<c)| (1<<d)| (1<<e)| (1<<f)| (1<<g)| (1<<h);
                              HammingBallMask.push_back(mask);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>11){
        //radius 12
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        for(int f = e+1; f < codelength; f++){
                          for(int g = f+1; g < codelength; g++){
                            for(int h = g+1; h < codelength; h++){
                              for(int l = h+1; h < codelength; l++){
                                unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b)| (1<<c)| (1<<d)| (1<<e)| (1<<f)| (1<<g)| (1<<h)| (1<<l);
                                HammingBallMask.push_back(mask);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

      if(radius>12){
        //radius 13
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        for(int f = e+1; f < codelength; f++){
                          for(int g = f+1; g < codelength; g++){
                            for(int h = g+1; h < codelength; h++){
                              for(int l = h+1; h < codelength; l++){
                                for(int m = l+1; m < codelength; m++){
                                  unsigned int mask = (1<<i) | (1<<j) | (1<<k)| (1<<a)| (1<<b)| (1<<c)| (1<<d)| (1<<e)| (1<<f)| (1<<g)| (1<<h)| (1<<l)| (1<<m);
                                  HammingBallMask.push_back(mask);
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask.size());
      }

    }

    void BuildHashTable64(int upbits, int lowbits, std::vector<Codes64>& baseAll ,std::vector<HashTable64>& tbAll){

      for(size_t h=0; h < baseAll.size(); h++){
        Codes64& base = baseAll[h];

        HashTable64 tb;
        for(int i = 0; i < (1 << upbits); i++){
          HashBucket64 emptyBucket;
          tb.push_back(emptyBucket);
        }

        for(size_t i = 0; i < base.size(); i ++){
          unsigned int idx1 = base[i] >> lowbits;
          unsigned long idx2 = base[i] - ((unsigned long)idx1 << lowbits);
          if(tb[idx1].find(idx2) != tb[idx1].end()){
            tb[idx1][idx2].push_back(i);
          }else{
            std::vector<unsigned int> v;
            v.push_back(i);
            tb[idx1].insert(make_pair(idx2,v));
          }
        }
        tbAll.push_back(tb);
      }
    }



    void generateMask64(){
      //i = 0 means the origin code
      HammingBallMask64.push_back(0);
      HammingRadius.push_back(HammingBallMask64.size());

      unsigned long One = 1;
      if(radius>0){
        //radius 1
        for(int i = 0; i < codelength; i++){
          unsigned long mask = One << i;
          HammingBallMask64.push_back(mask);
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>1){
        //radius 2
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            unsigned long mask = (One<<i) | (One<<j);
            HammingBallMask64.push_back(mask);
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>2){
        //radius 3
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              unsigned long mask = (One<<i) | (One<<j) | (One<<k);
              HammingBallMask64.push_back(mask);
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>3){
        //radius 4
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                unsigned long mask = (One<<i) | (One<<j) | (One<<k)| (One<<a);
                HammingBallMask64.push_back(mask);
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>4){
        //radius 5
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  unsigned long mask = (One<<i) | (One<<j) | (One<<k)| (One<<a)| (One<<b);
                  HammingBallMask64.push_back(mask);
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>5){
        //radius 6
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    unsigned long mask = (One<<i) | (One<<j) | (One<<k)| (One<<a)| (One<<b)| (One<<c);
                    HammingBallMask64.push_back(mask);
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>6){
        //radius 7
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      unsigned long mask = (One<<i) | (One<<j) | (One<<k)| (One<<a)| (One<<b)| (One<<c)| (One<<d);
                      HammingBallMask64.push_back(mask);
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>7){
        //radius 8
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        unsigned long mask = (One<<i) | (One<<j) | (One<<k)| (One<<a)| (One<<b)| (One<<c)| (One<<d)| (One<<e);
                        HammingBallMask64.push_back(mask);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>8){
        //radius 9
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        for(int f = e+1; f < codelength; f++){
                          unsigned long mask = (One<<i) | (One<<j) | (One<<k)| (One<<a)| (One<<b)| (One<<c)| (One<<d)| (One<<e)| (One<<f);
                          HammingBallMask64.push_back(mask);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>9){
        //radius 10
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        for(int f = e+1; f < codelength; f++){
                          for(int g = f+1; g < codelength; g++){
                            unsigned long mask = (One<<i) | (One<<j) | (One<<k)| (One<<a)| (One<<b)| (One<<c)| (One<<d)| (One<<e)| (One<<f)| (One<<g);
                            HammingBallMask64.push_back(mask);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

      if(radius>10){
        //radius 11
        for(int i = 0; i < codelength; i++){
          for(int j = i+1; j < codelength; j++){
            for(int k = j+1; k < codelength; k++){
              for(int a = k+1; a < codelength; a++){
                for(int b = a+1; b < codelength; b++){
                  for(int c = b+1; c < codelength; c++){
                    for(int d = c+1; d < codelength; d++){
                      for(int e = d+1; e < codelength; e++){
                        for(int f = e+1; f < codelength; f++){
                          for(int g = f+1; g < codelength; g++){
                            for(int h = g+1; h < codelength; h++){
                              unsigned long mask = (One<<i) | (One<<j) | (One<<k)| (One<<a)| (One<<b)| (One<<c)| (One<<d)| (One<<e)| (One<<f)| (One<<g)| (One<<h);
                              HammingBallMask64.push_back(mask);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        HammingRadius.push_back(HammingBallMask64.size());
      }

    }

    void buildIndexImpl()
    {
      std::cout<<"HASHING building hashing table"<<std::endl;

      if (codelength <= 32 ){
        codelength = codelength - codelengthshift;
        BuildHashTable32(upbits, codelength-upbits, BaseCode ,htb);
        generateMask32();
      }else if(codelength <= 64 ){
        codelength = codelength - codelengthshift;
        BuildHashTable64(upbits, codelength-upbits, BaseCode64 ,htb64);
        generateMask64();
      }else{
        std::cout<<"code length not supported yet!"<<std::endl;
      }
    }

    void getNeighbors(size_t K, const Matrix<DataType>& query){

      if(gs.size() != features_.get_rows()){
        if (codelength <= 32 ){
          getNeighbors32(K,query);
        }else if(codelength <= 64 ){
          getNeighbors64(K,query);
        }else{
          std::cout<<"code length not supported yet!"<<std::endl;
        }
      }else{
        switch(SP.search_method){
        case 0:
          if (codelength <= 32 ){
            getNeighborsIEH32_kgraph(K, query);
          }else if(codelength <= 64 ){
            getNeighborsIEH64_kgraph(K, query);
          }else{
            std::cout<<"code length not supported yet!"<<std::endl;
          }
          break;
        case 1:
          if (codelength <= 32 ){
            getNeighborsIEH32_nnexp(K, query);
          }else if(codelength <= 64 ){
            getNeighborsIEH64_nnexp(K, query);
          }else{
            std::cout<<"code length not supported yet!"<<std::endl;
          }
          break;
        default:
          std::cout<<"no such searching method"<<std::endl;
        }
      }

    }


    void getNeighbors32(size_t K, const Matrix<DataType>& query){
      int lowbits = codelength - upbits;

      unsigned int MaxCheck=HammingRadius[radius];
      std::cout<<"maxcheck : "<<MaxCheck<<std::endl;

      boost::dynamic_bitset<> tbflag(features_.get_rows(), false);

      nn_results.clear();

      VisitBucketNum.clear();
      VisitBucketNum.resize(radius+2);

      for(size_t cur = 0; cur < query.get_rows(); cur++){

        std::vector<unsigned int> pool(SP.search_init_num);
        unsigned int p = 0;
        tbflag.reset();

        unsigned int j = 0;
        for(; j < MaxCheck; j++){
          for(unsigned int h=0; h < QueryCode.size(); h++){
            unsigned int searchcode = QueryCode[h][cur] ^ HammingBallMask[j];
            unsigned int idx1 = searchcode >> lowbits;
            unsigned int idx2 = searchcode - (idx1 << lowbits);

            HashBucket::iterator bucket= htb[h][idx1].find(idx2);
            if(bucket != htb[h][idx1].end()){
              std::vector<unsigned int> vp = bucket->second;
              for(size_t k = 0; k < vp.size() && p < (unsigned int)SP.search_init_num; k++){
                if(tbflag.test(vp[k]))continue;

                tbflag.set(vp[k]);
                pool[p++]=(vp[k]);
              }
              if(p >= (unsigned int)SP.search_init_num)  break;
            }
            if(p >= (unsigned int)SP.search_init_num)  break;
          }
          if(p >= (unsigned int)SP.search_init_num)  break;
        }

        if(p < (unsigned int)SP.search_init_num){
          VisitBucketNum[radius+1]++;
        }else{
          for(int r=0;r<=radius;r++){
            if(j<=HammingRadius[r]){
              VisitBucketNum[r]++;
              break;
            }
          }
        }

        if (p<K){
          int base_n = features_.get_rows();
          while(p < K){
            unsigned int nn = rand() % base_n;
            if(tbflag.test(nn)) continue;
            tbflag.set(nn);
            pool[p++] = (nn);
          }

        }

        std::vector<std::pair<float,size_t>> result;
        for(unsigned int i=0; i<p;i++){
          result.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(pool[i]), features_.get_cols()),pool[i]));
        }
        std::partial_sort(result.begin(), result.begin() + K, result.end());

        std::vector<int> res;
        for(unsigned int j = 0; j < K; j++) res.push_back(result[j].second);
        nn_results.push_back(res);

      }
      //std::cout<<"bad query number:  " << VisitBucketNum[radius+1] << std::endl;
    }

    void getNeighbors64(size_t K, const Matrix<DataType>& query){
      int lowbits = codelength - upbits;

      unsigned int MaxCheck=HammingRadius[radius];
      std::cout<<"maxcheck : "<<MaxCheck<<std::endl;

      boost::dynamic_bitset<> tbflag(features_.get_rows(), false);

      nn_results.clear();

      VisitBucketNum.clear();
      VisitBucketNum.resize(radius+2);


      for(size_t cur = 0; cur < query.get_rows(); cur++){

        std::vector<unsigned int> pool(SP.search_init_num);
        unsigned int p = 0;
        tbflag.reset();

        unsigned int j = 0;
        for(; j < MaxCheck; j++){
          for(unsigned int h=0; h < QueryCode64.size(); h++){
            unsigned long searchcode = QueryCode64[h][cur] ^ HammingBallMask64[j];
            unsigned int idx1 = searchcode >> lowbits;
            unsigned long idx2 = searchcode - (( unsigned long)idx1 << lowbits);

            HashBucket64::iterator bucket= htb64[h][idx1].find(idx2);
            if(bucket != htb64[h][idx1].end()){
              std::vector<unsigned int> vp = bucket->second;
              for(size_t k = 0; k < vp.size() && p < (unsigned int)SP.search_init_num; k++){
                if(tbflag.test(vp[k]))continue;

                tbflag.set(vp[k]);
                pool[p++]=(vp[k]);
              }
              if(p >= (unsigned int)SP.search_init_num)  break;
            }
            if(p >= (unsigned int)SP.search_init_num)  break;
          }
          if(p >= (unsigned int)SP.search_init_num)  break;
        }


        if(p < (unsigned int)SP.search_init_num){
          VisitBucketNum[radius+1]++;
        }else{
          for(int r=0;r<=radius;r++){
            if(j<=HammingRadius[r]){
              VisitBucketNum[r]++;
              break;
            }
          }
        }


        if (p<K){
          int base_n = features_.get_rows();
          while(p < K){
            unsigned int nn = rand() % base_n;
            if(tbflag.test(nn)) continue;
            tbflag.set(nn);
            pool[p++] = (nn);
          }

        }

        std::vector<std::pair<float,size_t>> result;
        for(unsigned int i=0; i<p;i++){
          result.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(pool[i]), features_.get_cols()),pool[i]));
        }
        std::partial_sort(result.begin(), result.begin() + K, result.end());

        std::vector<int> res;
        for(unsigned int j = 0; j < K; j++) res.push_back(result[j].second);
        nn_results.push_back(res);
      }
      //std::cout<<"bad query number:  " <<VisitBucketNum[radius+1]<< std::endl;
    }

    void getNeighborsIEH32_nnexp(size_t K, const Matrix<DataType>& query){
      int lowbits = codelength - upbits;

      unsigned int MaxCheck=HammingRadius[radius];
      std::cout<<"maxcheck : "<<MaxCheck<<std::endl;

      int resultSize = SP.extend_to;
      if (K > (unsigned)SP.extend_to)
        resultSize = K;

      boost::dynamic_bitset<> tbflag(features_.get_rows(), false);
      nn_results.clear();

      VisitBucketNum.clear();
      VisitBucketNum.resize(radius+2);

      for(size_t cur = 0; cur < query.get_rows(); cur++){
        tbflag.reset();

        std::vector<int> pool(SP.search_init_num);
        unsigned int p = 0;


        unsigned int j = 0;
        for(; j < MaxCheck; j++){
          for(size_t h=0; h < QueryCode.size(); h++){
            unsigned int searchcode = QueryCode[h][cur] ^ HammingBallMask[j];
            unsigned int idx1 = searchcode >> lowbits;
            unsigned int idx2 = searchcode - (idx1 << lowbits);

            HashBucket::iterator bucket= htb[h][idx1].find(idx2);
            if(bucket != htb[h][idx1].end()){
              std::vector<unsigned int> vp = bucket->second;
              for(size_t k = 0; k < vp.size() && p < (unsigned int)SP.search_init_num; k++){
                if(tbflag.test(vp[k]))continue;

                tbflag.set(vp[k]);
                pool[p++]=(vp[k]);
              }
              if(p >= (unsigned int)SP.search_init_num)  break;
            }
            if(p >= (unsigned int)SP.search_init_num)  break;
          }
          if(p >= (unsigned int)SP.search_init_num)  break;
        }

        if(p < (unsigned int)SP.search_init_num){
          VisitBucketNum[radius+1]++;
        }else{
          for(int r=0;r<=radius;r++){
            if(j<=HammingRadius[r]){
              VisitBucketNum[r]++;
              break;
            }
          }
        }

        int base_n = features_.get_rows();
        while(p < (unsigned int)SP.search_init_num){
          unsigned int nn = rand() % base_n;
          if(tbflag.test(nn)) continue;
          tbflag.set(nn);
          pool[p++] = (nn);
        }

        //sorting the pool
        std::vector<std::pair<float,size_t>> result;
        for(unsigned int i=0; i<pool.size();i++){
          result.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(pool[i]), features_.get_cols()),pool[i]));
        }
        std::partial_sort(result.begin(), result.begin() + resultSize, result.end());
        result.resize(resultSize);
        pool.clear();
        for(int j = 0; j < resultSize; j++)
          pool.push_back(result[j].second);

        //nn_exp
        boost::dynamic_bitset<> newflag(features_.get_rows(), true);
        newflag.set();

        int iter=0;
        std::vector<int> ids;
        while(iter++ < SP.search_epoches){
          //the heap is max heap
          ids.clear();
          for(unsigned j = 0; j < SP.extend_to ; j++){
            if(newflag.test( pool[j] )){
              newflag.reset(pool[j]);

              for(unsigned neighbor=0; neighbor < gs[pool[j]].size(); neighbor++){
                unsigned id = gs[pool[j]][neighbor];

                if(tbflag.test(id))continue;
                else tbflag.set(id);

                ids.push_back(id);
              }
            }
          }

          for(size_t j = 0; j < ids.size(); j++){
            result.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(ids[j]), features_.get_cols()),ids[j]));
          }
          std::partial_sort(result.begin(), result.begin() + resultSize, result.end());
          result.resize(resultSize);
          pool.clear();
          for(int j = 0; j < resultSize; j++)
            pool.push_back(result[j].second);
        }

        if(K<(unsigned)SP.extend_to)
          pool.resize(K);

        nn_results.push_back(pool);
      }
    }

    void getNeighborsIEH32_kgraph(size_t K, const Matrix<DataType>& query){
      int lowbits = codelength - upbits;

      unsigned int MaxCheck=HammingRadius[radius];
      std::cout<<"maxcheck : "<<MaxCheck<<std::endl;

      nn_results.clear();
      boost::dynamic_bitset<> tbflag(features_.get_rows(), false);

      bool bSorted = true;
      unsigned pool_size = SP.search_epoches * SP.extend_to;
      if (pool_size >= (unsigned)SP.search_init_num){
        SP.search_init_num = pool_size;
        bSorted = false;
      }

      VisitBucketNum.clear();
      VisitBucketNum.resize(radius+2);


      for(size_t cur = 0; cur < query.get_rows(); cur++){

        tbflag.reset();
        std::vector<unsigned int> pool(SP.search_init_num);
        unsigned int p = 0;

        unsigned int j = 0;
        for(; j < MaxCheck; j++){
          for(size_t h=0; h < QueryCode.size(); h++){
            unsigned int searchcode = QueryCode[h][cur] ^ HammingBallMask[j];
            unsigned int idx1 = searchcode >> lowbits;
            unsigned int idx2 = searchcode - (idx1 << lowbits);

            HashBucket::iterator bucket= htb[h][idx1].find(idx2);
            if(bucket != htb[h][idx1].end()){
              std::vector<unsigned int> vp = bucket->second;
              for(size_t k = 0; k < vp.size() && p < (unsigned int)SP.search_init_num; k++){
                if(tbflag.test(vp[k]))continue;

                tbflag.set(vp[k]);
                pool[p++]=(vp[k]);
              }
              if(p >= (unsigned int)SP.search_init_num)  break;
            }
            if(p >= (unsigned int)SP.search_init_num)  break;
          }
          if(p >= (unsigned int)SP.search_init_num)  break;
        }

        if(p < (unsigned int)SP.search_init_num){
          VisitBucketNum[radius+1]++;
        }else{
          for(int r=0;r<=radius;r++){
            if(j<=HammingRadius[r]){
              VisitBucketNum[r]++;
              break;
            }
          }
        }

        int base_n = features_.get_rows();
        while(p < (unsigned int)SP.search_init_num){
          unsigned int nn = rand() % base_n;
          if(tbflag.test(nn)) continue;
          tbflag.set(nn);
          pool[p++] = (nn);
        }

        std::vector<std::pair<float,size_t>> result;
        for(unsigned int i=0; i<pool.size();i++){
          result.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(pool[i]), features_.get_cols()),pool[i]));
        }
        if(bSorted){
          std::partial_sort(result.begin(), result.begin() + pool_size, result.end());
          result.resize(pool_size);
        }


        tbflag.reset();

        std::vector<Point> knn(K + SP.extend_to +1);
        std::vector<Point> results;

        for (unsigned iter = 0; iter < (unsigned)SP.search_epoches; iter++) {

          unsigned L = 0;
          for(unsigned j=0; j < SP.extend_to ; j++){
            if(!tbflag.test(result[iter*SP.extend_to+j].second)){
              tbflag.set(result[iter*SP.extend_to+j].second);
              knn[L].id = result[iter*SP.extend_to+j].second;
              knn[L].dist = result[iter*SP.extend_to+j].first;
              knn[L].flag = true;
              L++;
            }
          }
          if(~bSorted){
            std::sort(knn.begin(), knn.begin() + L);
          }

          unsigned int k =  0;
          while (k < L) {
            unsigned int nk = L;
            if (knn[k].flag) {
              knn[k].flag = false;
              unsigned n = knn[k].id;

              for(unsigned neighbor=0; neighbor < gs[n].size(); neighbor++){
                unsigned id = gs[n][neighbor];

                if(tbflag.test(id))continue;
                tbflag.set(id);

                float dist = distance_->compare(query.get_row(cur), features_.get_row(id), features_.get_cols());
                Point nn(id, dist);
                unsigned int r = InsertIntoKnn(&knn[0], L, nn);

                //if ( (r <= L) && (L + 1 < knn.size())) ++L;
                if ( L + 1 < knn.size()) ++L;
                if (r < nk) nk = r;
              }
            }
            if (nk <= k) k = nk;
            else ++k;
          }


          if (L > K) L = K;
          if (results.empty()) {
            results.reserve(K + 1);
            results.resize(L + 1);
            std::copy(knn.begin(), knn.begin() + L, results.begin());
          }
          else {
            for (unsigned int l = 0; l < L; ++l) {
              unsigned r = InsertIntoKnn(&results[0], results.size() - 1, knn[l]);
              if (r < results.size() /* inserted */ && results.size() < (K + 1)) {
                results.resize(results.size() + 1);
              }
            }
          }
        }

        std::vector<int> res;
        for(size_t i = 0; i < K && i < results.size();i++)
          res.push_back(results[i].id);
        nn_results.push_back(res);

      }
    }

    void getNeighborsIEH64_nnexp(size_t K, const Matrix<DataType>& query){

      int lowbits = codelength - upbits;

      unsigned int MaxCheck=HammingRadius[radius];
      std::cout<<"maxcheck : "<<MaxCheck<<std::endl;

      int resultSize = SP.extend_to;
      if (K > (unsigned)SP.extend_to)
        resultSize = K;

      boost::dynamic_bitset<> tbflag(features_.get_rows(), false);
      nn_results.clear();

      VisitBucketNum.clear();
      VisitBucketNum.resize(radius+2);

      for(size_t cur = 0; cur < query.get_rows(); cur++){

        std::vector<int> pool(SP.search_init_num);
        unsigned int p = 0;
        tbflag.reset();

        unsigned int j = 0;
        for(; j < MaxCheck; j++){
          for(unsigned int h=0; h < QueryCode64.size(); h++){
            unsigned long searchcode = QueryCode64[h][cur] ^ HammingBallMask64[j];
            unsigned int idx1 = searchcode >> lowbits;
            unsigned long idx2 = searchcode - (( unsigned long)idx1 << lowbits);

            HashBucket64::iterator bucket= htb64[h][idx1].find(idx2);
            if(bucket != htb64[h][idx1].end()){
              std::vector<unsigned int> vp = bucket->second;
              for(size_t k = 0; k < vp.size() && p < (unsigned int)SP.search_init_num; k++){
                if(tbflag.test(vp[k]))continue;

                tbflag.set(vp[k]);
                pool[p++]=(vp[k]);
              }
              if(p >= (unsigned int)SP.search_init_num)  break;
            }
            if(p >= (unsigned int)SP.search_init_num)  break;
          }
          if(p >= (unsigned int)SP.search_init_num)  break;
        }

        if(p < (unsigned int)SP.search_init_num){
          VisitBucketNum[radius+1]++;
        }else{
          for(int r=0;r<=radius;r++){
            if(j<=HammingRadius[r]){
              VisitBucketNum[r]++;
              break;
            }
          }
        }

        int base_n = features_.get_rows();
        while(p < (unsigned int)SP.search_init_num){
          unsigned int nn = rand() % base_n;
          if(tbflag.test(nn)) continue;
          tbflag.set(nn);
          pool[p++] = (nn);
        }

        //sorting the pool
        std::vector<std::pair<float,size_t>> result;
        for(unsigned int i=0; i<pool.size();i++){
          result.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(pool[i]), features_.get_cols()),pool[i]));
        }
        std::partial_sort(result.begin(), result.begin() + resultSize, result.end());
        result.resize(resultSize);
        pool.clear();
        for(int j = 0; j < resultSize; j++)
          pool.push_back(result[j].second);

        //nn_exp
        boost::dynamic_bitset<> newflag(features_.get_rows(), true);
        newflag.set();

        int iter=0;
        std::vector<int> ids;
        while(iter++ < SP.search_epoches){
          //the heap is max heap
          ids.clear();
          for(unsigned j = 0; j < SP.extend_to ; j++){
            if(newflag.test( pool[j] )){
              newflag.reset(pool[j]);

              for(unsigned neighbor=0; neighbor < gs[pool[j]].size(); neighbor++){
                unsigned id = gs[pool[j]][neighbor];

                if(tbflag.test(id))continue;
                else tbflag.set(id);

                ids.push_back(id);
              }
            }
          }

          for(size_t j = 0; j < ids.size(); j++){
            result.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(ids[j]), features_.get_cols()),ids[j]));
          }
          std::partial_sort(result.begin(), result.begin() + resultSize, result.end());
          result.resize(resultSize);
          pool.clear();
          for(int j = 0; j < resultSize; j++)
            pool.push_back(result[j].second);
        }

        if(K<(unsigned)SP.extend_to)
          pool.resize(K);

        nn_results.push_back(pool);
      }
    }

    void getNeighborsIEH64_kgraph(size_t K, const Matrix<DataType>& query){
      int lowbits = codelength - upbits;

      unsigned int MaxCheck=HammingRadius[radius];
      std::cout<<"maxcheck : "<<MaxCheck<<std::endl;

      nn_results.clear();
      boost::dynamic_bitset<> tbflag(features_.get_rows(), false);

      bool bSorted = true;
      unsigned pool_size = SP.search_epoches * SP.extend_to;
      if (pool_size >= (unsigned)SP.search_init_num){
        SP.search_init_num = pool_size;
        bSorted = false;
      }

      VisitBucketNum.clear();
      VisitBucketNum.resize(radius+2);


      for(size_t cur = 0; cur < query.get_rows(); cur++){

        std::vector<unsigned int> pool(SP.search_init_num);
        unsigned int p = 0;
        tbflag.reset();

        unsigned int j = 0;
        for(; j < MaxCheck; j++){
          for(unsigned int h=0; h < QueryCode64.size(); h++){
            unsigned long searchcode = QueryCode64[h][cur] ^ HammingBallMask64[j];
            unsigned int idx1 = searchcode >> lowbits;
            unsigned long idx2 = searchcode - (( unsigned long)idx1 << lowbits);

            HashBucket64::iterator bucket= htb64[h][idx1].find(idx2);
            if(bucket != htb64[h][idx1].end()){
              std::vector<unsigned int> vp = bucket->second;
              for(size_t k = 0; k < vp.size() && p < (unsigned int)SP.search_init_num; k++){
                if(tbflag.test(vp[k]))continue;

                tbflag.set(vp[k]);
                pool[p++]=(vp[k]);
              }
              if(p >= (unsigned int)SP.search_init_num)  break;
            }
            if(p >= (unsigned int)SP.search_init_num)  break;
          }
          if(p >= (unsigned int)SP.search_init_num)  break;
        }

        if(p < (unsigned int)SP.search_init_num){
          VisitBucketNum[radius+1]++;
        }else{
          for(int r=0;r<=radius;r++){
            if(j<=HammingRadius[r]){
              VisitBucketNum[r]++;
              break;
            }
          }
        }

        int base_n = features_.get_rows();
        while(p < (unsigned int)SP.search_init_num){
          unsigned int nn = rand() % base_n;
          if(tbflag.test(nn)) continue;
          tbflag.set(nn);
          pool[p++] = (nn);
        }

        std::vector<std::pair<float,size_t>> result;
        for(unsigned int i=0; i<pool.size();i++){
          result.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(pool[i]), features_.get_cols()),pool[i]));
        }
        if(bSorted){
          std::partial_sort(result.begin(), result.begin() + pool_size, result.end());
          result.resize(pool_size);
        }


        tbflag.reset();

        std::vector<Point> knn(K + SP.extend_to +1);
        std::vector<Point> results;

        for (unsigned iter = 0; iter < (unsigned)SP.search_epoches; iter++) {

          unsigned L = 0;
          for(unsigned j=0; j < (unsigned)SP.extend_to ; j++){
            if(!tbflag.test(result[iter*SP.extend_to+j].second)){
              tbflag.set(result[iter*SP.extend_to+j].second);
              knn[L].id = result[iter*SP.extend_to+j].second;
              knn[L].dist = result[iter*SP.extend_to+j].first;
              knn[L].flag = true;
              L++;
            }
          }
          if(~bSorted){
            std::sort(knn.begin(), knn.begin() + L);
          }

          unsigned int k =  0;
          while (k < L) {
            unsigned int nk = L;
            if (knn[k].flag) {
              knn[k].flag = false;
              unsigned n = knn[k].id;

              for(unsigned neighbor=0; neighbor < gs[n].size(); neighbor++){
                unsigned id = gs[n][neighbor];

                if(tbflag.test(id))continue;
                tbflag.set(id);

                float dist = distance_->compare(query.get_row(cur), features_.get_row(id), features_.get_cols());
                Point nn(id, dist);
                unsigned int r = InsertIntoKnn(&knn[0], L, nn);

                //if ( (r <= L) && (L + 1 < knn.size())) ++L;
                if ( L + 1 < knn.size()) ++L;
                if (r < nk) nk = r;
              }
            }
            if (nk <= k) k = nk;
            else ++k;
          }


          if (L > K) L = K;
          if (results.empty()) {
            results.reserve(K + 1);
            results.resize(L + 1);
            std::copy(knn.begin(), knn.begin() + L, results.begin());
          }
          else {
            for (unsigned int l = 0; l < L; ++l) {
              unsigned r = InsertIntoKnn(&results[0], results.size() - 1, knn[l]);
              if (r < results.size() /* inserted */ && results.size() < (K + 1)) {
                results.resize(results.size() + 1);
              }
            }
          }
        }

        std::vector<int> res;
        for(size_t i = 0; i < K && i < results.size();i++)
          res.push_back(results[i].id);
        nn_results.push_back(res);

      }
    }

    void outputVisitBucketNum(){
      unsigned i=0;
      std::cout<< "Radius " << i <<" bucket num: "<<HammingRadius[i]<< " points num: "<< VisitBucketNum[i]<<std::endl;

      for(i=1; i<HammingRadius.size();i++){
        std::cout<< "Radius " << i <<" bucket num: "<<HammingRadius[i] - HammingRadius[i-1]<< " points num: "<< VisitBucketNum[i]<<std::endl;
      }
      std::cout<< "Radius larger, points num: " << VisitBucketNum[i]<<std::endl;
    }

    void loadIndex(char* filename){}
    void saveIndex(char* filename){}
    void loadTrees(char* filename){}
    void saveTrees(char* filename){}
    void loadGraph(char* filename){
      std::ifstream in(filename,std::ios::binary);
      in.seekg(0,std::ios::end);
      std::ios::pos_type ss = in.tellg();
      size_t fsize = (size_t)ss;
      int dim;
      in.seekg(0,std::ios::beg);
      in.read((char*)&dim, sizeof(int));
      size_t num = fsize / (dim+1) / 4;
      std::cout<<"load g "<<num<<" "<<dim<< std::endl;
      in.seekg(0,std::ios::beg);
      gs.clear();
      for(size_t i = 0; i < num; i++){
        std::vector<unsigned> heap;
        in.read((char*)&dim, sizeof(int));
        for(int j =0; j < dim; j++){
          unsigned id;
          in.read((char*)&id, sizeof(int));
          heap.push_back(id);
        }
        gs.push_back(heap);
      }
      in.close();
    }
    void saveGraph(char* filename){}
    void initGraph(){}

  protected:
    int tablenum;
    int upbits;
    int codelength;
    int codelengthshift;
    int radius;
    USING_BASECLASS_SYMBOLS

    std::vector<HashTable> htb;
    std::vector<Codes> BaseCode;
    std::vector<Codes> QueryCode;
    std::vector<unsigned int> HammingBallMask;

    std::vector<HashTable64> htb64;
    std::vector<Codes64> BaseCode64;
    std::vector<Codes64> QueryCode64;
    std::vector<unsigned long> HammingBallMask64;

    std::vector<unsigned int> HammingRadius;

    // for statistic info
    std::vector<unsigned int> VisitBucketNum;

};
}
#endif

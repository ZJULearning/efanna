#ifndef EFANNA_IEH_INDEX_H_
#define EFANNA_IEH_INDEX_H_

#include "algorithm/base_index.hpp"
#include <boost/dynamic_bitset.hpp>
#include <time.h>
//for Debug
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#define MAX_ROWSIZE 1024
#define HASH_RADIUS 1
#define DEPTH 16
#define MAX_POOL_SIZE 5500


namespace efanna{
  struct IEHIndexParams : public IndexParams
  {
      IEHIndexParams(int UpperBits, int HashRadius, char*& BaseCodeFile, char*& QueryCodeFile)
      {
          init_index_type = IEH;
          ValueType bcf;
          bcf.str_pt = BaseCodeFile;
          extra_params.insert(std::make_pair("bcfile",bcf));
          ValueType qcf;
          qcf.str_pt = QueryCodeFile;
          extra_params.insert(std::make_pair("qcfile",qcf));
          ValueType upb;
          upb.int_val = UpperBits;
          extra_params.insert(std::make_pair("upbits",upb));
          ValueType radius;
          radius.int_val = HashRadius;
          extra_params.insert(std::make_pair("radius",radius));
      }
  };

  template <typename DataType>
  class IEHIndex : public InitIndex<DataType>
  {
  public:

    typedef InitIndex<DataType> BaseClass;
    typedef std::vector<unsigned int> Codes;
    typedef std::unordered_map<unsigned int, std::vector<unsigned int> > HashBucket;
    typedef std::vector<HashBucket> HashTable;
    IEHIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = IEHIndexParams(0,NULL,NULL)) :
      BaseClass(dataset,d,params)
    {
      std::cout<<"IEH initial, max code length : 32, max hasing radius : 6."<<std::endl;
      ExtraParamsMap::const_iterator it = params_.extra_params.find("bcfile");
      if(it != params_.extra_params.end()){
        char* fpath = (it->second).str_pt;
        std::string str(fpath);
        std::cout << "Loading base code from " << str << std::endl;
        LoadCode(fpath, BaseCode);
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
        LoadCode(fpath, QueryCode);
        std::cout << "code length is "<<codelength<<std::endl;
      }
      else{
        std::cout << "error: no query code file" << std::endl;
      }

      it = params_.extra_params.find("upbits");
      upbits = (it->second).int_val;
      if(it != params_.extra_params.end()){
        std::cout << "use upper "<<upbits<< " bits as first level index of hashtable"<< std::endl;
        std::cout << "use lower "<<codelength - upbits<< " bits as second level index of hashtable"<< std::endl;
      }
      else{
        std::cout << "error: no upper bits number setting" << std::endl;
      }

      it = params_.extra_params.find("radius");
      radius = (it->second).int_val;
      if(it != params_.extra_params.end()){
        std::cout << "search hamming radius "<<radius<< std::endl;
        if(radius > 6){
          std::cout << "radius greater than 6 not supported yet!" << std::endl;
          radius = 6;
        }
      }
      else{
        std::cout << "error: no radius number setting" << std::endl;
      }
    }



    /**
     * Builds the index
     */
     void StringSplit(std::string src, std::vector<std::string>& des){
       int start = 0;
       int end = 0;
       for(size_t i = 0; i < src.length(); i++){
         if(src[i]==' '){
           end = i;
           //if(end>start)cout<<start<<" "<<end<<" "<<src.substr(start,end-start)<<endl;
           des.push_back(src.substr(start,end-start));
           start = i+1;
         }
       }
     }

     void LoadCode(char* filename, Codes& base){
       std::ifstream in(filename);
       char buf[MAX_ROWSIZE];
       //int cnt = 0;
       while(!in.eof()){
         in.getline(buf,MAX_ROWSIZE);
         std::string strtmp(buf);
         std::vector<std::string> strs;
         StringSplit(strtmp,strs);
         if(strs.size()<2)continue;
         unsigned int codetmp = 0;
         codelength = strs.size();
         for(size_t i = 0; i < strs.size(); i++){
           unsigned int c = atoi(strs[i].c_str());
           codetmp = codetmp << 1;
           codetmp += c;

         }//if(cnt++ > 999998){cout<<strs.size()<<" "<<buf<<" "<<codetmp<<endl;}
         base.push_back(codetmp);
       }std::cout<<"point num: "<<base.size()<<std::endl;
       in.close();
     }

     void BuildHashTable(int upbits, int lowbits, Codes base ,HashTable& tb){
       tb.clear();
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
     }

     void generateMask(){
       //i = 0 means the origin code
       HammingBallMask.push_back(0);
       //radius 1
       for(int i = 0; i < codelength; i++){
         unsigned int mask = 1 << i;
         HammingBallMask.push_back(mask);
       }
       //radius 2
       for(int i = 0; i < codelength; i++){
         for(int j = i+1; j < codelength; j++){
           unsigned int mask = (1<<i) | (1<<j);
           HammingBallMask.push_back(mask);
         }
       }
       //radius 3
       for(int i = 0; i < codelength; i++){
         for(int j = i+1; j < codelength; j++){
             for(int k = j+1; k < codelength; k++){
                 unsigned int mask = (1<<i) | (1<<j) | (1<<k);
                 HammingBallMask.push_back(mask);
             }
         }
       }
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
        //std::cout << HammingBallMask.size() << std::endl;
       //for(size_t i = 0; i < HammingBallMask.size(); i++)std::cout << HammingBallMask[i]<< " ";
     }

     void buildIndexImpl()
     {
       std::cout<<"IEH building hashing table"<<std::endl;
       BuildHashTable(upbits, codelength - upbits, BaseCode ,htb);
       generateMask();
     }


     void getNeighbors_version1(size_t k, const Matrix<DataType>& query){
       int lowbits = codelength - upbits;
       int MaxCheck=1;

       switch (radius) {
         case 0:
           MaxCheck = 1;break;
         case 1:
           MaxCheck = codelength * 1 + 1;break;
         case 2:
           MaxCheck = codelength * (codelength + 1)/2 + 1;break;
         case 3:
           MaxCheck = codelength * (codelength - 1)* (codelength - 2)/6 + codelength * (codelength - 1)/2 + codelength * 1 + 1;break;
       }


       boost::dynamic_bitset<> tbflag(features_.get_rows(), false);
       boost::dynamic_bitset<> newflag(features_.get_rows(), false);

       std::cout<<"maxcheck : "<<MaxCheck<<std::endl;

       nn_results.clear();


       for(size_t cur = 0; cur < query.get_rows(); cur++){

         CandidateHeap Candidates;
         std::vector<int> candstmp;
         tbflag.reset();
         newflag.reset();

         for(int j = 0; j < MaxCheck; j++){
           unsigned int searchcode = QueryCode[cur] ^ HammingBallMask[j];
           unsigned int idx1 = searchcode >> lowbits;
           unsigned int idx2 = searchcode - (idx1 << lowbits);

           HashBucket::iterator bucket= htb[idx1].find(idx2);
           if(bucket != htb[idx1].end()){
             std::vector<unsigned int> vp = bucket->second;
             for(size_t k = 0; k < vp.size() && candstmp.size() < (unsigned int)SP.search_init_num; k++){
               candstmp.push_back(vp[k]);
             }
             if(candstmp.size() >= (unsigned int)SP.search_init_num){
              break;
             }
           }
           if(candstmp.size() >= (unsigned int)SP.search_init_num){
            break;
           }
         }
//std::cout<<candstmp.size()<<std::endl;

         for(size_t j = 0; j < candstmp.size(); j++){
           int nn = candstmp[j];

           newflag.set(nn);
           tbflag.set(nn);

           Candidate<DataType> c(nn, distance_->compare(
                 query.get_row(cur), features_.get_row(nn), features_.get_cols())
               );
           Candidates.insert(c);
           if(Candidates.size() > (unsigned int)SP.search_init_num)
              Candidates.erase(Candidates.begin());
         }



         int base_n = features_.get_rows();

   while(Candidates.size() < (unsigned int)SP.search_init_num){
            unsigned int nn = rand() % base_n;
            if(tbflag.test(nn))continue;
            else{ tbflag.set(nn); newflag.set(nn);}

            Candidate<DataType> c(nn, distance_->compare(
                  query.get_row(cur), features_.get_row(nn), features_.get_cols())
                );
          Candidates.insert(c);
         }

     int iter=0;
         std::vector<int> ids;
         while(iter++ < SP.search_epoches){
           //the heap is max heap
           typename CandidateHeap::reverse_iterator it = Candidates.rbegin();
           ids.clear();
           for(int j = 0; j < SP.extend_to && it != Candidates.rend(); j++,it++){
             if(newflag.test(it->row_id)){
               newflag.reset(it->row_id);

               typename CandidateHeap::reverse_iterator neighbor = knn_graph[it->row_id].rbegin();

               for(size_t nnk = 0;nnk < knn_graph[it->row_id].size() && neighbor != knn_graph[it->row_id].rend(); neighbor++, nnk++){
                 if(tbflag.test(neighbor->row_id))continue;
                 else tbflag.set(neighbor->row_id);

                  ids.push_back(neighbor->row_id);
               }
             }
           }
           for(size_t j = 0; j < ids.size(); j++){

             Candidate<DataType> c(ids[j], distance_->compare(
                 query.get_row(cur), features_.get_row(ids[j]), features_.get_cols())
               );

             Candidates.insert(c);
             newflag.set(ids[j]);
             if(Candidates.size() > (unsigned int)SP.search_init_num)Candidates.erase(Candidates.begin());
           }
         }
         typename CandidateHeap::reverse_iterator it = Candidates.rbegin();
         std::vector<int> res;
         for(size_t i = 0; i < k && it != Candidates.rend();i++, it++)
            res.push_back(it->row_id);
         nn_results.push_back(res);
       }

     }

     void getNeighbors(size_t K, const Matrix<DataType>& query){
       int lowbits = codelength - upbits;
       int MaxCheck=1;
       if (0<radius)  MaxCheck = MaxCheck +  codelength;
       if (1<radius)  MaxCheck = MaxCheck +  codelength * (codelength - 1)/2;
       if (2<radius)  MaxCheck = MaxCheck +   codelength * (codelength - 1)* (codelength - 2)/6 ;
       if (3<radius)  MaxCheck = MaxCheck +   codelength * (codelength - 1)* (codelength - 2)* (codelength - 3)/24;
       if (4<radius)  MaxCheck = MaxCheck +  codelength * (codelength - 1)* (codelength - 2)* (codelength - 3)* (codelength - 4)/120 ;
       if (5<radius)  MaxCheck = MaxCheck +  codelength * (codelength - 1)* (codelength - 2)* (codelength - 3)* (codelength - 4)* (codelength - 4)/720 ;
       if (6<radius)  return;

       std::cout<<"maxcheck : "<<MaxCheck<<std::endl;


       boost::dynamic_bitset<> tbflag(features_.get_rows(), false);
       nn_results.clear();

       for(size_t cur = 0; cur < query.get_rows(); cur++){
         tbflag.reset();

         std::vector<unsigned int> pool(SP.search_init_num);
         unsigned int p = 0;

         for(int j = 0; j < MaxCheck; j++){
           unsigned int searchcode = QueryCode[cur] ^ HammingBallMask[j];
           unsigned int idx1 = searchcode >> lowbits;
           unsigned int idx2 = searchcode - (idx1 << lowbits);

           HashBucket::iterator bucket= htb[idx1].find(idx2);
           if(bucket != htb[idx1].end()){
             std::vector<unsigned int> vp = bucket->second;
             for(size_t k = 0; k < vp.size() && p < (unsigned int)SP.search_init_num; k++){
            	 pool[p++]=(vp[k]);
            	 tbflag.set(vp[k]);
             }
             if(p >= (unsigned int)SP.search_init_num)  break;
           }
           if(p >= (unsigned int)SP.search_init_num)  break;
         }

         int base_n = features_.get_rows();
         while(p < (unsigned int)SP.search_init_num){
           unsigned int nn = rand() % base_n;
           if(tbflag.test(nn)) continue;
           pool[p++] = (nn);
         }

         std::vector<int> res;
         InitIndex<DataType>::nnExpansion(K, query.get_row(cur), pool, res);
         nn_results.push_back(res);
       }
     }


  void getNeighbors_kgraph(size_t K, const Matrix<DataType>& query){
       int lowbits = codelength - upbits;
       int MaxCheck=1;

       switch (radius) {
         case 0:
           MaxCheck = 1;break;
         case 1:
           MaxCheck = codelength * 1 + 1;break;
         case 2:
           MaxCheck = codelength * (codelength + 1)/2 + 1;break;
         case 3:
           MaxCheck = codelength * (codelength - 1)* (codelength - 2)/6 + codelength * (codelength - 1)/2 + codelength * 1 + 1;break;
       }

       std::cout<<"maxcheck : "<<MaxCheck<<std::endl;
       //clock_t sum = 0;

       nn_results.clear();
       boost::dynamic_bitset<> flag(features_.get_rows(), false);

       for(size_t cur = 0; cur < query.get_rows(); cur++){

         flag.reset();
         std::vector<unsigned int> pool(SP.search_init_num);
         unsigned int p = 0;

         for(int j = 0; j < MaxCheck; j++){
           unsigned int searchcode = QueryCode[cur] ^ HammingBallMask[j];
           unsigned int idx1 = searchcode >> lowbits;
           unsigned int idx2 = searchcode - (idx1 << lowbits);

           HashBucket::iterator bucket= htb[idx1].find(idx2);
           if(bucket != htb[idx1].end()){
             std::vector<unsigned int> vp = bucket->second;
             for(size_t k = 0; k < vp.size() && p < (unsigned int)SP.search_init_num; k++){
            	 pool[p++]=(vp[k]);
               flag.set(vp[k]);
             }
             if(p >= (unsigned int)SP.search_init_num){
              break;
             }
           }
           if(p >= (unsigned int)SP.search_init_num){
            break;
           }
         }
         //std::cout<<candstmp.size()<<std::endl;

         int base_n = features_.get_rows();

         while(p < (unsigned int)SP.search_init_num){
           unsigned int nn = rand() % base_n;
           if(flag.test(nn)) continue;
           pool[p++] = (nn);
         }

         std::vector<Point> results;

         InitIndex<DataType>::nnExpansion_kgraph(K, query.get_row(cur), pool, results);

         std::vector<int> res;
         for(size_t i = 0; i < K && i < results.size();i++)
            res.push_back(results[i].id);
         nn_results.push_back(res);
       }
       //std::cout<<"Call nnExpansion time : "<<sum/CLOCKS_PER_SEC<<" seconds"<<std::endl;
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
       knn_graph.clear();
       for(size_t i = 0; i < num; i++){
         CandidateHeap heap;
         in.read((char*)&dim, sizeof(int));
         for(int j =0; j < dim; j++){
           int id;
           in.read((char*)&id, sizeof(int));
           Candidate<DataType> can(id, -1);
           heap.insert(can);
         }
         knn_graph.push_back(heap);
       }
       in.close();
     }
     void saveGraph(char* filename){}
     void initGraph(){}
     typedef std::set<Candidate<DataType>, std::greater<Candidate<DataType>> > CandidateHeap;

   protected:
     int upbits;
     int codelength;
     int radius;
     HashTable htb;
     USING_BASECLASS_SYMBOLS
     Codes BaseCode;
     Codes QueryCode;
     std::vector<unsigned int> HammingBallMask;

   };
}
#endif

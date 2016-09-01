#ifndef EFANNA_NNEXP_INDEX_H_
#define EFANNA_NNEXP_INDEX_H_

#include "algorithm/base_index.hpp"
#include <algorithm>


#include <time.h>



//for Debug
#include <iostream>
namespace efanna{
  struct NNexpIndexParams : public IndexParams
  {
      NNexpIndexParams(int epoches = 10, int extend = 10, int building_use_k = 10)
      {
          init_index_type = NNEXP;
          K = building_use_k;
          build_epoches = epoches;
          extend_num = extend;
      }
  };

  template <typename DataType>
  class NNexpIndex : public InitIndex<DataType>
  {
  public:

    typedef InitIndex<DataType> BaseClass;
    NNexpIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = NNexpIndexParams()) :
      BaseClass(dataset,d,params)
    {
      std::cout<<"NN-expansion initial"<<std::endl;
    }



    /**
     * Builds the index
     */

     void buildIndexImpl()
     {
       std::cout<<"NN-expansion build"<<std::endl;
       initGraph();
       int iter = 0;
       while(iter < params_.build_epoches){
         NNexpBuild();iter++;
         std::cout<<"iteration "<<iter<<std::endl;
       }
     }

     //void SetZeros(std::vector<unsigned char>& v){
     //    for(size_t i=0; i < v.size(); i++)v[i]=0;
     //}


     void getNeighbors_version1(size_t k, const Matrix<DataType>& query){

       int base_n = features_.get_rows();
       boost::dynamic_bitset<> tbflag(base_n, false);
       boost::dynamic_bitset<> newflag(base_n, false);

       nn_results.clear();


       for(unsigned int cur = 0; cur < query.get_rows(); cur++){

         CandidateHeap Candidates;
         tbflag.reset();
         newflag.reset();

         while(Candidates.size() < (unsigned int)SP.search_init_num){
           unsigned int nn = rand() % base_n;
           Candidate<DataType> c(nn, distance_->compare(
                 query.get_row(cur), features_.get_row(nn), features_.get_cols())
               );
           Candidates.insert(c);
           newflag.set(nn);
           tbflag.set(nn);
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
                for(size_t nnk = 0;nnk < params_.K && neighbor != knn_graph[it->row_id].rend(); neighbor++, nnk++){
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

         int base_n = features_.get_rows();
         nn_results.clear();


         for(unsigned int cur = 0; cur < query.get_rows(); cur++){

             std::vector<unsigned int> pool(SP.search_init_num);
             for(unsigned int i = 0;i< (unsigned int) SP.search_init_num;i++){
            	 pool[i] = rand() % base_n;
             }

             std::vector<int> res;
             nnExpansion(K, query.get_row(cur), pool, res);
             nn_results.push_back(res);

         }

    }

     void getNeighbors_kgraph(size_t K, const Matrix<DataType>& query){

         int base_n = features_.get_rows();
         nn_results.clear();

         //clock_t s,f,sum = 0;
         //clock_t sum = 0;

         for(unsigned int cur = 0; cur < query.get_rows(); cur++){

             std::vector<unsigned int> pool(SP.search_init_num);
             for(unsigned int i = 0;i< (unsigned int) SP.search_init_num;i++){
            	 pool[i] = rand() % base_n;
             }

            // s = clock();

             std::vector<Point> results;
             nnExpansion_kgraph(K,  query.get_row(cur), pool, results);

            // f = clock();
             //sum = sum + f-s;

             std::vector<int> res;
             for(size_t i = 0; i < K && i < results.size();i++)
                res.push_back(results[i].id);
             nn_results.push_back(res);
         }
        //std::cout<<"Call nnExpansion time : "<<sum/CLOCKS_PER_SEC<<" seconds"<<std::endl;

    }


     void loadIndex(char* filename){}
     void saveIndex(char* filename){
        std::ofstream out(filename,std::ios::binary);
      /*std::vector<std::vector<int>>::iterator i;//naive implemented
       for(i = knn_graph.begin(); i!= knn_graph.end(); i++){
         std::vector<int>::iterator j;
         int dim = params_.K;
         out.write((char*)&dim, sizeof(int));
         for(j = i->begin(); j != i->end(); j++){
           int id = *j;
           out.write((char*)&id, sizeof(int));
         }
       }*/

       int dim = params_.K;
       for(size_t i = 0; i < knn_graph.size(); i++){
         typename CandidateHeap::reverse_iterator it = knn_graph[i].rbegin();
         out.write((char*)&dim, sizeof(int));
         for(size_t j=0 ; j<params_.K && it!= knn_graph[i].rend(); j++,it++ ){
           int id = it->row_id;
           out.write((char*)&id, sizeof(int));
         }
       }
       out.close();
     }
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
      //std::cout<<"load g "<<num<<" "<<dim<< std::endl;
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
    void saveGraph(char* filename){
     std::ofstream out(filename,std::ios::binary);

     int dim = params_.K;
     for(size_t i = 0; i < knn_graph.size(); i++){
       typename CandidateHeap::reverse_iterator it = knn_graph[i].rbegin();
       out.write((char*)&dim, sizeof(int));
       for(size_t j =0; j < params_.K && it!= knn_graph[i].rend(); j++,it++ ){
         int id = it->row_id;
         out.write((char*)&id, sizeof(int));
       }
     }

     out.close();
    }
    static void GenRandom (unsigned *addr, unsigned size, unsigned N) {
        for (unsigned i = 0; i < size; ++i) {
            addr[i] = rand() % (N - size);
        }
        std::sort(addr, addr + size);
        for (unsigned i = 1; i < size; ++i) {
            if (addr[i] <= addr[i-1]) {
                addr[i] = addr[i-1] + 1;
            }
        }
        unsigned off = rand() % N;
        for (unsigned i = 0; i < size; ++i) {
            addr[i] = (addr[i] + off) % N;
        }
    }
    void initGraph(){
      int N = features_.get_rows();

      std::vector<unsigned> random(params_.extend_num);
//random initialize
      for (unsigned n = 0; n < (unsigned)N; ++n) {
          GenRandom(&random[0], random.size(), N);
          CandidateHeap can;
          for(size_t j = 0; j < random.size(); j++){
            DataType dist = distance_->compare(
              features_.get_row(n), features_.get_row(random[j]), features_.get_cols());
            Candidate<DataType> c(random[j],dist);
            can.insert(c);
          }
          knn_graph.push_back(can);
      }
    }


    void NNexpBuild(){
      //int tenp = knn_graph.size() / 10;
      size_t base_n = knn_graph.size();
      int *check = new int[base_n];
      size_t y;


      for(size_t x = 0; x < base_n; ++x) check[x] = -1;
      for(size_t i = 0; i < base_n; i++){
        //if(i % tenp == 0)std::cout<<i / tenp<<" / 10 completed"<<std::endl;
        typename CandidateHeap::reverse_iterator rit = knn_graph[i].rbegin();
        check[i] = i;

        std::vector<size_t> ids;
        ids.clear();
        rit = knn_graph[i].rbegin();
        for(;rit != knn_graph[i].rend(); rit++){
          typename CandidateHeap::reverse_iterator neighbor = knn_graph[rit->row_id].rbegin();
          y = rit->row_id;
          check[y] = i;
          for(; neighbor != knn_graph[rit->row_id].rend(); neighbor++){
            ids.push_back(neighbor->row_id);
          }
        }
        for(size_t j = 0; j < ids.size(); j++)if(check[ y=ids[j] ]!=(int)i){
          check[y] = i;
          DataType dist =  distance_->compare(features_.get_row(i), features_.get_row(ids[j]), features_.get_cols());
          if(dist < knn_graph[i].begin()->distance){
            Candidate<DataType> c(ids[j], dist );
            knn_graph[i].insert(c);
            if(knn_graph[i].size() > (unsigned)params_.extend_num)knn_graph[i].erase(knn_graph[i].begin());
          }

          if(dist < knn_graph[ids[j]].begin()->distance){
            Candidate<DataType> c(i, dist );
            knn_graph[ids[j]].insert(c);
            if(knn_graph[ids[j]].size() > (unsigned int)params_.extend_num)knn_graph[ids[j]].erase(knn_graph[ids[j]].begin());
          }

        }

      }delete[] check;
    }
    typedef std::set<Candidate<DataType>, std::greater<Candidate<DataType>> > CandidateHeap;
  protected:
    USING_BASECLASS_SYMBOLS
  };
}
#endif

#ifndef EFANNA_KDTREE_UB_INDEX_H_
#define EFANNA_KDTREE_UB_INDEX_H_
#include "algorithm/base_index.hpp"
#include <fstream>
#include <ctime>
#include <string.h>
#include <random>
#include <queue>
//#include <bitset>
//using std::bitset;
#include <boost/dynamic_bitset.hpp>

namespace efanna{
struct KDTreeUbIndexParams : public IndexParams
{
	KDTreeUbIndexParams(bool rnn_used, int tree_num_total, int merge_level = 4, int epoches = 4, int check = 25, int myL = 30, int building_use_k = 10, int tree_num_build = 0, int myS = 10)
	{
		reverse_nn_used = rnn_used;
		init_index_type = KDTREE_UB;
		K = building_use_k;
		build_epoches = epoches;
		S = myS;
		ValueType treev;
		treev.int_val = tree_num_total;
		extra_params.insert(std::make_pair("trees",treev));
		ValueType treeb;
		treeb.int_val = tree_num_build > 0 ? tree_num_build : tree_num_total;
		extra_params.insert(std::make_pair("treesb",treeb));
		ValueType merge_levelv;
		merge_levelv.int_val = merge_level;
		extra_params.insert(std::make_pair("ml",merge_levelv));
		L = myL;
		Check_K = check;
	}
};
template <typename DataType>
class KDTreeUbIndex : public InitIndex<DataType>
{
public:

	typedef InitIndex<DataType> BaseClass;
	KDTreeUbIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = KDTreeUbIndexParams(true,4)) :
		BaseClass(dataset,d,params)
	{
		std::cout<<"kdtree ub initial"<<std::endl;
		ExtraParamsMap::const_iterator it = params_.extra_params.find("trees");
		if(it != params_.extra_params.end()){
			TreeNum = (it->second).int_val;
#ifdef INFO
			std::cout << "Using kdtree to build "<< TreeNum << " trees in total" << std::endl;
#endif
		}
		else{
			TreeNum = 4;
#ifdef INFO
			std::cout << "Using kdtree to build "<< TreeNum << " trees in total" << std::endl;
#endif
		}
		SP.tree_num = TreeNum;

		it = params_.extra_params.find("treesb");
		if(it != params_.extra_params.end()){
			TreeNumBuild = (it->second).int_val;
#ifdef INFO
			std::cout << "Building kdtree graph with "<< TreeNumBuild <<" trees"<< std::endl;
#endif
		}
		else{
			TreeNumBuild = TreeNum;
#ifdef INFO
			std::cout << "Building kdtree graph with "<< TreeNumBuild <<" trees"<< std::endl;
#endif
		}

		it = params_.extra_params.find("ml");
		if(it != params_.extra_params.end()){
			ml = (it->second).int_val;
#ifdef INFO
			std::cout << "Building kdtree initial index with merge level "<< ml  << std::endl;
#endif
		}
		else{
			ml = -1;
#ifdef INFO
			std::cout << "Building kdtree initial index with max merge level "<< std::endl;
#endif
		}
		max_deepth = 0x0fffffff;
		error_flag = false;
	}

	void buildIndexImpl(){
#ifdef INFO
		clock_t s,f;
		s = clock();
#endif
		initGraph();

#ifdef INFO
		f = clock();
#endif

		std::cout << "initial graph finised"<< std::endl;
#ifdef INFO
		std::cout << "initial graph using time: "<< (f-s)*1.0/CLOCKS_PER_SEC<<" seconds"<< std::endl;
#endif

		if(error_flag){
			std::cout << "merge level deeper than tree, max merge deepth is" << max_deepth-1<<std::endl;
			return;
		}
		refineGraph();
	}
	struct Node
	{
		int DivDim;
		DataType DivVal;
		size_t StartIdx, EndIdx;
		unsigned treeid;
		Node* Lchild, * Rchild;

		~Node() {
			if (Lchild!=NULL) Lchild->~Node();
			if (Rchild!=NULL) Rchild->~Node();
		}

	};

	void loadIndex(char* filename){
		read_data(filename);
	}
	void saveIndex(char* filename){

		size_t points_num = features_.get_rows();
		size_t feature_dim = features_.get_cols();
		save_data(filename, params_.K, points_num, feature_dim);
	}
	//algorithms copy and rewrite from flann
	void loadTrees(char* filename){
		std::ifstream in(filename, std::ios::binary|std::ios::in);
		if(!in.is_open()){std::cout<<"open file error"<<std::endl;exit(-10087);}
		unsigned int K,tree_num;
		size_t dim,num;

		//read file head
		in.read((char*)&(K),sizeof(unsigned int));
		in.read((char*)&(tree_num),sizeof(unsigned int));
		in.read((char*)&(num),sizeof(size_t));
		in.read((char*)&(dim),sizeof(size_t));

		SP.tree_num = tree_num;

		//read trees

		tree_roots_.clear();
		for(unsigned int i=0;i<tree_num;i++){// for each tree
			int node_num, node_size;
			in.read((char*)&(node_num),sizeof(int));
			in.read((char*)&(node_size),sizeof(int));

			std::vector<struct Node *> tree_nodes;
			for(int j=0;j<node_num;j++){
				struct Node *tmp = new struct Node();
				in.read((char*)&(tmp->DivDim),sizeof(tmp->DivDim));
				in.read((char*)&(tmp->DivVal),sizeof(tmp->DivVal));
				in.read((char*)&(tmp->StartIdx),sizeof(tmp->StartIdx));
				in.read((char*)&(tmp->EndIdx),sizeof(tmp->EndIdx));
				in.read((char*)&(tmp->Lchild),sizeof(tmp->Lchild));
				in.read((char*)&(tmp->Rchild),sizeof(tmp->Rchild));
				tmp->Lchild = NULL;
				tmp->Rchild = NULL;
				tmp->treeid = i;
				tree_nodes.push_back(tmp);


			}
			//std::cout<<"build "<<i<<std::endl;
			struct Node *root = DepthFirstBuildTree(tree_nodes);
			if(root==NULL){ exit(-11); }
			tree_roots_.push_back(root);
		}

		//read index range
		LeafLists.clear();
		for(unsigned int i=0;i<tree_num;i++){

			std::vector<unsigned> leaves;
			for(unsigned int j=0;j<num; j++){
				unsigned leaf;
				in.read((char*)&(leaf),sizeof(int));
				leaves.push_back(leaf);
			}
			LeafLists.push_back(leaves);
		}
		in.close();
	}
	void saveTrees(char* filename){
		unsigned int K = params_.K;
		size_t num = features_.get_rows();
		size_t dim = features_.get_cols();
		std::fstream out(filename, std::ios::binary|std::ios::out);
		if(!out.is_open()){std::cout<<"open file error"<<std::endl;exit(-10086);}
		unsigned int tree_num = tree_roots_.size();

		//write file head
		out.write((char *)&(K), sizeof(unsigned int));
		out.write((char *)&(tree_num), sizeof(unsigned int));
		out.write((char *)&(num), sizeof(size_t)); //feature point number
		out.write((char *)&(dim), sizeof(size_t)); //feature dim

		//write trees
		typename std::vector<Node *>::iterator it;//int cnt=0;
		for(it=tree_roots_.begin(); it!=tree_roots_.end(); it++){
			//write tree nodes with depth first trace


			size_t offset_node_num = out.tellp();

			out.seekp(sizeof(int),std::ios::cur);

			unsigned int node_size = sizeof(struct Node);
			out.write((char *)&(node_size), sizeof(int));

			unsigned int node_num = DepthFirstWrite(out, *it);

			out.seekg(offset_node_num,std::ios::beg);

			out.write((char *)&(node_num), sizeof(int));

			out.seekp(0,std::ios::end);
			//std::cout<<"tree: "<<cnt++<<" written, node: "<<node_num<<" at offset " << offset_node_num <<std::endl;
		}

		if(LeafLists.size()!=tree_num){ std::cout << "leaf_size!=tree_num" << std::endl; exit(-6); }

		for(unsigned int i=0; i<tree_num; i++){
			for(unsigned int j=0;j<num;j++){
				out.write((char *)&(LeafLists[i][j]), sizeof(int));
			}
		}
		out.close();
	}
	void loadGraph(char* filename){
		std::ifstream in(filename,std::ios::binary);
		unsigned N;

		in.seekg(0,std::ios::end);
		std::ios::pos_type ss = in.tellg();
		size_t fsize = (size_t)ss;
		int dim;
		in.seekg(0,std::ios::beg);
		in.read((char*)&dim, sizeof(int));
		N = fsize / (dim+1) / 4;

		in.seekg(0,std::ios::beg);

		gs.resize(N);
		//M.resize(N);
		//norms.resize(N);
		for(unsigned i=0; i < N; i++){
			unsigned k;
			//DataType norm;
			in.read((char*)&k, sizeof(unsigned));
			//in.read((char*)&m, sizeof(unsigned));
			//in.read((char*)&norm, sizeof(DataType));
			//norms[i] = norm;
			//M[i] = m;
			gs[i].resize(k);

			for(unsigned j=0; j<k; j++){
				unsigned id;
				in.read((char*)&id, sizeof(unsigned));
				gs[i][j] = id;
			}
		}
		in.close();
	}
	/*
    void saveGraph(char* filename){
     std::ofstream out(filename,std::ios::binary);

     int dim = params_.K;//int meansize = 0;
     for(size_t i = 0; i < knn_graph.size(); i++){
       typename CandidateHeap::reverse_iterator it = knn_graph[i].rbegin();
       out.write((char*)&dim, sizeof(int));//meansize += knn_graph[i].size();
       for(size_t j =0; j < params_.K && it!= knn_graph[i].rend(); j++,it++ ){
         int id = it->row_id;
         out.write((char*)&id, sizeof(int));
       }
     }//meansize /= knn_graph.size();
     //std::cout << "size mean " << meansize << std::endl;
     out.close();
    }
	 */
	void saveGraph(char* filename){
		std::ofstream out(filename,std::ios::binary);
		unsigned N = gs.size();
		//out.write((char*)&N, sizeof(int));
		for(unsigned i=0; i < N; i++){
			unsigned k = gs[i].size();
			//unsigned m = M[i];
			//DataType norm = norms[i];
			out.write((char*)&k, sizeof(unsigned));
			//out.write((char*)&m, sizeof(unsigned));
			//out.write((char*)&norm, sizeof(DataType));
			for(unsigned j = 0; j < k; j++){
				unsigned id = gs[i][j];
				out.write((char*)&id, sizeof(unsigned));
			}
		}
		out.close();
	}
	//for nn search

	void SearchQueryToLeaf(Node* node, const DataType* q, unsigned dep, std::vector<Node*>& node_pool){
		if(node->Lchild != NULL && node->Rchild !=NULL){
			if(q[node->DivDim] < node->DivVal){
				SearchQueryToLeaf(node->Lchild, q, dep, node_pool);
				if(node_pool.size() < dep)
					SearchQueryToLeaf(node->Rchild, q, dep, node_pool);
			}
			else{
				SearchQueryToLeaf(node->Rchild, q, dep, node_pool);
				if(node_pool.size() < dep)
					SearchQueryToLeaf(node->Lchild, q, dep, node_pool);
			}
		}
		else
			node_pool.push_back(node);
	}

	void getSearchNodeList(Node* node, const DataType* q, unsigned int lsize, std::vector<Node*>& vn){
		if(vn.size() >= lsize)
			return;

		if(node->Lchild != NULL && node->Rchild !=NULL){
			if(q[node->DivDim] < node->DivVal){
				getSearchNodeList(node->Lchild, q, lsize,  vn );
				getSearchNodeList(node->Rchild, q, lsize, vn);
			}else{
				getSearchNodeList(node->Rchild, q, lsize, vn);
				getSearchNodeList(node->Lchild, q, lsize, vn);
			}
		}else
			vn.push_back(node);
	}


	void getNeighbors(size_t searchK, const Matrix<DataType>& query){
		switch(SP.search_method){
		case 0:
			getNeighbors_nnexp(searchK, query);
			break;
		case 1:
			getNeighbors_kgraph(searchK, query);
			break;
		default:
			std::cout<<"no such searching method"<<std::endl;
		}

	}

	void getNeighbors_nnexp(size_t K, const Matrix<DataType>& query){
#ifdef INFO
		std::cout<<"using tree num "<< SP.tree_num<<std::endl;
#endif
		if(SP.tree_num > tree_roots_.size()){
			std::cout<<"wrong tree number"<<std::endl;return;
		}

		nn_results.clear();
		nn_results.resize(query.get_rows());
		unsigned dim = features_.get_cols();

		int resultSize = SP.extend_to;
		if (K > (unsigned)SP.extend_to)
			resultSize = K;


#pragma omp parallel for
		for(unsigned int cur = 0; cur < query.get_rows(); cur++){
			boost::dynamic_bitset<> tbflag(features_.get_rows(), false);
			boost::dynamic_bitset<> newflag(features_.get_rows(), true);
			tbflag.reset();
			newflag.set();

			std::vector<std::vector<Node*>> NodeCandi;
			NodeCandi.resize(SP.tree_num);

			const DataType* q_row = query.get_row(cur);
			_mm_prefetch((char *)q_row, _MM_HINT_T0);
			unsigned int lsize = SP.search_init_num*2 / (5*SP.tree_num) + 1;
			for(unsigned int i = 0; i < SP.tree_num; i++){
				getSearchNodeList(tree_roots_[i], q_row, lsize, NodeCandi[i]);
			}
			std::vector<int> pool(SP.search_init_num);
			unsigned int p = 0;
			for(unsigned int ni = 0; ni < lsize; ni++){
				for(unsigned int i = 0; i < NodeCandi.size(); i++){
					Node* leafn = NodeCandi[i][ni];
					for(size_t j = leafn->StartIdx; j < leafn->EndIdx && p < (unsigned int)SP.search_init_num; j++){
						size_t nn = LeafLists[i][j];
						if(tbflag.test(nn))continue;
						tbflag.set(nn);
						pool[p++]=(nn);
					}
					if(p >= (unsigned int)SP.search_init_num) break;
				}
				if(p >= (unsigned int)SP.search_init_num) break;
			}
			int base_n = features_.get_rows();
			while(p < (unsigned int)SP.search_init_num){
				unsigned int nn = rand() % base_n;
				if(tbflag.test(nn))continue;
				tbflag.set(nn);
				pool[p++]=(nn);
			}


			std::vector<std::pair<float,size_t>> result;
			//for(unsigned int i=0; i<pool.size();i++){
			//  _mm_prefetch((char *)features_.get_row(pool[i]), _MM_HINT_T0);
			//}
			unsigned cache_blocksz = 80;
			for(unsigned int i=0; i*cache_blocksz<pool.size();i++){
				unsigned s = i*cache_blocksz;
				unsigned t = s + cache_blocksz > pool.size() ? pool.size() : s+cache_blocksz;
				unsigned s_ = s;
				while(s<t){
					_mm_prefetch((char *)features_.get_row(pool[s]), _MM_HINT_T0);
					s++;
				}
				while(s_<t){
					result.push_back(std::make_pair(distance_->compare(q_row, features_.get_row(pool[s_]), dim),pool[s_]));
					s_++;
				}
			}
			std::partial_sort(result.begin(), result.begin() + resultSize, result.end());
			result.resize(resultSize);
			pool.clear();
			for(int j = 0; j < resultSize; j++)
				pool.push_back(result[j].second);

			int iter=0;
			std::vector<int> ids;
			while(iter++ < SP.search_epoches){
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
				//for(unsigned int j=0; j<ids.size();j++){
				//_mm_prefetch((char *)features_.get_row(ids[j]), _MM_HINT_T0);
				//}
				for(size_t j = 0; j * cache_blocksz< ids.size(); j++){
					unsigned s = j * cache_blocksz;
					unsigned t = s + cache_blocksz > ids.size() ? ids.size() : s+cache_blocksz;
					unsigned s_ = s;
					while(s<t){
						_mm_prefetch((char *)features_.get_row(ids[s]), _MM_HINT_T0);
						s++;
					}
					while(s_<t){
						result.push_back(std::make_pair(distance_->compare(q_row, features_.get_row(ids[s_]), dim),ids[s_]));
						s_++;
					}
					//result.push_back(std::make_pair(distance_->compare(q_row, features_.get_row(ids[j]), dim),ids[j]));
				}
				std::partial_sort(result.begin(), result.begin() + resultSize, result.end());
				result.resize(resultSize);
				pool.clear();
				for(int j = 0; j < resultSize; j++)
					pool.push_back(result[j].second);
			}

			if(K<SP.extend_to)
				pool.resize(K);

			//nn_results.push_back(pool);
			std::vector<int>& res = nn_results[cur];
			for(unsigned i = 0; i < K ;i++)
				res.push_back(pool[i]);
		}
	}

	void getNeighbors_kgraph(size_t searchK, const Matrix<DataType>& query){
#ifdef INFO
		std::cout<<"using tree num "<< SP.tree_num<<std::endl;
#endif
		if(SP.tree_num > tree_roots_.size()){
			std::cout<<"wrong tree number"<<std::endl;return;
		}

		nn_results.clear();
		nn_results.resize(query.get_rows());
		unsigned dim = features_.get_cols();
		unsigned int lsize = SP.search_init_num*2 / (5*SP.tree_num) + 1;

		bool bSorted = true;
		unsigned pool_size = SP.search_epoches * SP.extend_to;
		if (pool_size >= (unsigned)SP.search_init_num){
			SP.search_init_num = pool_size;
			bSorted = false;
		}

#pragma omp parallel for
		for(unsigned int cur = 0; cur < query.get_rows(); cur++){
			std::mt19937 rng(1998);
			boost::dynamic_bitset<> flags(features_.get_rows(), false);

			std::vector<std::vector<Node*> > Vnl;
			Vnl.resize(SP.tree_num);
			const DataType* q_row = query.get_row(cur);
			_mm_prefetch((char *)q_row, _MM_HINT_T0);
			for(unsigned int i = 0; i < SP.tree_num; i++){
				getSearchNodeList(tree_roots_[i], q_row, lsize, Vnl[i]);
			}

			std::vector<int> pool(SP.search_init_num);
			unsigned int p = 0;
			for(unsigned int ni = 0; ni < lsize; ni++){
				for(unsigned int i = 0; i < Vnl.size(); i++){
					Node* leafn = Vnl[i][ni];
					for(size_t j = leafn->StartIdx; j < leafn->EndIdx && p < (unsigned int)SP.search_init_num; j++){
						size_t nn = LeafLists[i][j];
						if(flags.test(nn))continue;
						flags.set(nn);
						pool[p++]=(nn);
					}
					if(p >= (unsigned int)SP.search_init_num) break;
				}
				if(p >= (unsigned int)SP.search_init_num) break;
			}
			int base_n = features_.get_rows();
			while(p < (unsigned int)SP.search_init_num){
				unsigned int nn = rand() % base_n;
				if(flags.test(nn))continue;
				flags.set(nn);
				pool[p++]=(nn);
			}

			std::vector<std::pair<float,size_t>> result;
			unsigned cache_blocksz = 80;
			for(unsigned int i=0; i*cache_blocksz<pool.size();i++){
				unsigned s = i*cache_blocksz;
				unsigned t = s + cache_blocksz > pool.size() ? pool.size() : s+cache_blocksz;
				unsigned s_ = s;
				while(s<t){
					_mm_prefetch((char *)features_.get_row(pool[s]), _MM_HINT_T0);
					s++;
				}
				while(s_<t){
					result.push_back(std::make_pair(distance_->compare(q_row, features_.get_row(pool[s_]), dim),pool[s_]));
					s_++;
				}
			}
			if(bSorted){
				std::partial_sort(result.begin(), result.begin() + pool_size, result.end());
				result.resize(pool_size);
			}

			flags.reset();
			std::vector<Point> knn(searchK + SP.extend_to +1);
			std::vector<Point> results;
			for (unsigned iter = 0; iter < (unsigned)SP.search_epoches; ++iter) {

				unsigned L = 0;
				for(unsigned j=0; j < (unsigned)SP.extend_to ; j++){
					if(!flags.test(result[iter*SP.extend_to+j].second)){
						flags.set(result[iter*SP.extend_to+j].second);
						knn[L].id = result[iter*SP.extend_to+j].second;
						knn[L].dist = result[iter*SP.extend_to+j].first;
						knn[L].flag = true;
						L++;
					}
				}
				if(~bSorted){
					std::sort(knn.begin(), knn.begin() + L);
				}

				unsigned k =  0;
				while (k < L) {
					unsigned nk = L;
					if (knn[k].flag) {
						knn[k].flag = false;
						unsigned n = knn[k].id;

						//unsigned maxM = M[n];
						unsigned maxM = SP.extend_to;
						//if ((unsigned)SP.extend_to > maxM) maxM = SP.extend_to;
						auto const &neighbors = gs[n];
						if (maxM > neighbors.size()) {
							maxM = neighbors.size();
						}

						for(unsigned m = 0; m < maxM; ++m){
							_mm_prefetch((char *)features_.get_row(neighbors[m]), _MM_HINT_T0);
						}
						for (unsigned m = 0; m < maxM; ++m) {
							unsigned id = neighbors[m];
							//BOOST_VERIFY(id < graph.size());
							if (flags[id]) continue;
							flags[id] = true;

							DataType dist = distance_->compare(q_row, features_.get_row(id), dim);

							Point nn(id, dist);
							unsigned r = InsertIntoKnn(&knn[0], L, nn);
							//BOOST_VERIFY(r <= L);
							//if (r > L) continue;
							if (L + 1 < knn.size()) ++L;
							if (r < nk) {
								nk = r;
							}
						}
					}
					if (nk <= k) {
						k = nk;
					}
					else {
						++k;
					}
				}
				if (L > searchK) L = searchK;

				if (results.empty()) {
					results.reserve(searchK + 1);
					results.resize(L + 1);
					std::copy(knn.begin(), knn.begin() + L, results.begin());
				} else {
					for (unsigned int l = 0; l < L; ++l) {
						unsigned r = InsertIntoKnn(&results[0], results.size() - 1, knn[l]);
						if (r < results.size()  && results.size() < (searchK + 1)) {
							results.resize(results.size() + 1);
						}
					}
				}
			}

			std::vector<int>& res = nn_results[cur];
			for(size_t i = 0; i < searchK && i < results.size();i++)
				res.push_back(results[i].id);
		}
	}



	int DepthFirstWrite(std::fstream& out, struct Node *root){
		if(root==NULL) return 0;
		int left_cnt = DepthFirstWrite(out, root->Lchild);
		int right_cnt = DepthFirstWrite(out, root->Rchild);

		//std::cout << root->StartIdx <<":" << root->EndIdx<< std::endl;
		out.write((char *)&(root->DivDim), sizeof(root->DivDim));
		out.write((char *)&(root->DivVal), sizeof(root->DivVal));
		out.write((char *)&(root->StartIdx), sizeof(root->StartIdx));
		out.write((char *)&(root->EndIdx), sizeof(root->EndIdx));
		out.write((char *)&(root->Lchild), sizeof(root->Lchild));
		out.write((char *)&(root->Rchild), sizeof(root->Rchild));
		return (left_cnt + right_cnt + 1);
	}
	struct Node* DepthFirstBuildTree(std::vector<struct Node *>& tree_nodes){
		std::vector<Node*> root_serial;
		typename std::vector<struct Node*>::iterator it = tree_nodes.begin();
		for( ; it!=tree_nodes.end(); it++){
			Node* tmp = *it;
			size_t rsize = root_serial.size();
			if(rsize<2){
				root_serial.push_back(tmp);
				//continue;
			}
			else{
				Node *last1 = root_serial[rsize-1];
				Node *last2 = root_serial[rsize-2];
				if(last1->EndIdx == tmp->EndIdx && last2->StartIdx == tmp->StartIdx){
					tmp->Rchild = last1;
					tmp->Lchild = last2;
					root_serial.pop_back();
					root_serial.pop_back();
				}
				root_serial.push_back(tmp);
			}

		}
		if(root_serial.size()!=1){
			std::cout << "Error constructing trees" << std::endl;
			return NULL;
		}
		return root_serial[0];
	}
	void read_data(char *filename){
		std::ifstream in(filename, std::ios::binary|std::ios::in);
		if(!in.is_open()){std::cout<<"open file error"<<std::endl;exit(-10087);}
		unsigned int K,tree_num;
		size_t dim,num;

		//read file head
		in.read((char*)&(K),sizeof(unsigned int));
		in.read((char*)&(tree_num),sizeof(unsigned int));
		in.read((char*)&(num),sizeof(size_t));
		in.read((char*)&(dim),sizeof(size_t));

		SP.tree_num = tree_num;

		//read trees

		tree_roots_.clear();
		for(unsigned int i=0;i<tree_num;i++){// for each tree
			int node_num, node_size;
			in.read((char*)&(node_num),sizeof(int));
			in.read((char*)&(node_size),sizeof(int));

			std::vector<struct Node *> tree_nodes;
			for(int j=0;j<node_num;j++){
				struct Node *tmp = new struct Node();
				in.read((char*)&(tmp->DivDim),sizeof(tmp->DivDim));
				in.read((char*)&(tmp->DivVal),sizeof(tmp->DivVal));
				in.read((char*)&(tmp->StartIdx),sizeof(tmp->StartIdx));
				in.read((char*)&(tmp->EndIdx),sizeof(tmp->EndIdx));
				in.read((char*)&(tmp->Lchild),sizeof(tmp->Lchild));
				in.read((char*)&(tmp->Rchild),sizeof(tmp->Rchild));
				tmp->Lchild = NULL;
				tmp->Rchild = NULL;
				tree_nodes.push_back(tmp);


			}
			//std::cout<<"build "<<i<<std::endl;
			struct Node *root = DepthFirstBuildTree(tree_nodes);
			if(root==NULL){ exit(-11); }
			tree_roots_.push_back(root);
		}

		//read index range
		LeafLists.clear();
		for(unsigned int i=0;i<tree_num;i++){

			std::vector<unsigned> leaves;
			for(unsigned int j=0;j<num; j++){
				unsigned leaf;
				in.read((char*)&(leaf),sizeof(int));
				leaves.push_back(leaf);
			}
			LeafLists.push_back(leaves);
		}

		//read knn graph
		knn_graph.clear();
		for(size_t i = 0; i < num; i++){
			CandidateHeap heap;
			for(size_t j =0; j < K ; j++ ){
				int id;
				in.read((char*)&id, sizeof(int));
				Candidate<DataType> can(id, -1);
				heap.insert(can);
			}
			knn_graph.push_back(heap);
		}
		in.close();
	}
	void save_data(char *filename, unsigned int K, size_t num, size_t dim){
		std::fstream out(filename, std::ios::binary|std::ios::out);
		if(!out.is_open()){std::cout<<"open file error"<<std::endl;exit(-10086);}
		unsigned int tree_num = tree_roots_.size();

		//write file head
		out.write((char *)&(K), sizeof(unsigned int));
		out.write((char *)&(tree_num), sizeof(unsigned int));
		out.write((char *)&(num), sizeof(size_t)); //feature point number
		out.write((char *)&(dim), sizeof(size_t)); //feature dim

		//write trees
		typename std::vector<Node *>::iterator it;//int cnt=0;
		for(it=tree_roots_.begin(); it!=tree_roots_.end(); it++){
			//write tree nodes with depth first trace


			size_t offset_node_num = out.tellp();

			out.seekp(sizeof(int),std::ios::cur);

			unsigned int node_size = sizeof(struct Node);
			out.write((char *)&(node_size), sizeof(int));

			unsigned int node_num = DepthFirstWrite(out, *it);

			out.seekg(offset_node_num,std::ios::beg);

			out.write((char *)&(node_num), sizeof(int));

			out.seekp(0,std::ios::end);
			//std::cout<<"tree: "<<cnt++<<" written, node: "<<node_num<<" at offset " << offset_node_num <<std::endl;
		}

		if(LeafLists.size()!=tree_num){ std::cout << "leaf_size!=tree_num" << std::endl; exit(-6); }

		for(unsigned int i=0; i<tree_num; i++){
			for(unsigned int j=0;j<num;j++){
				out.write((char *)&(LeafLists[i][j]), sizeof(int));
			}
		}

		//write knn-graph

		if(knn_graph.size()!=num){std::cout << "Error:" << std::endl; exit(-1);}
		for(size_t i = 0; i < knn_graph.size(); i++){
			typename CandidateHeap::reverse_iterator it = knn_graph[i].rbegin();
			for(size_t j =0; j < K && it!= knn_graph[i].rend(); j++,it++ ){
				int id = it->row_id;
				out.write((char*)&id, sizeof(int));
			}
		}

		out.close();
	}
	/*
    Node* divideTree(std::mt19937& rng, int* indices, size_t count, size_t offset){
      Node* node = new Node();
      if(count <= params_.TNS){
        node->DivDim = -1;
        node->Lchild = NULL;
        node->Rchild = NULL;
        node->StartIdx = offset;
        node->EndIdx = offset + count;
        //add points

        for(size_t i = 0; i < count; i++){
          for(size_t j = i+1; j < count; j++){
            DataType dist = distance_->compare(
                features_.get_row(indices[i]), features_.get_row(indices[j]), features_.get_cols());

            if(knn_graph[indices[i]].size() < params_.S || dist < knn_graph[indices[i]].begin()->distance){
              Candidate<DataType> c1(indices[j], dist);
              knn_graph[indices[i]].insert(c1);
              if(knn_graph[indices[i]].size() > params_.S)knn_graph[indices[i]].erase(knn_graph[indices[i]].begin());
            }
            else if(nhoods[indices[i]].nn_new.size() < params_.S * 2)nhoods[indices[i]].nn_new.push_back(indices[j]);
            if(knn_graph[indices[j]].size() < params_.S || dist < knn_graph[indices[j]].begin()->distance){
              Candidate<DataType> c2(indices[i], dist);
              knn_graph[indices[j]].insert(c2);
              if(knn_graph[indices[j]].size() > params_.S)knn_graph[indices[j]].erase(knn_graph[indices[j]].begin());
            }
            else if(nhoods[indices[j]].nn_new.size() < params_.S * 2)nhoods[indices[j]].nn_new.push_back(indices[i]);
          }
        }

      }else{
        int idx;
        int cutdim;
        DataType cutval;
        meanSplit(rng, indices, count, idx, cutdim, cutval);

        node->DivDim = cutdim;
        node->DivVal = cutval;
        node->StartIdx = offset;
        node->EndIdx = offset + count;
        node->Lchild = divideTree(rng, indices, idx, offset);
        node->Rchild = divideTree(rng, indices+idx, count-idx, offset+idx);
      }

      return node;
    }

    Node* divideTreeOnly(std::mt19937& rng, unsigned* indices, size_t count, size_t offset){
      Node* node = new Node();
      if(count <= params_.TNS){
        node->DivDim = -1;
        node->Lchild = NULL;
        node->Rchild = NULL;
        node->StartIdx = offset;
        node->EndIdx = offset + count;
        //add points

      }else{
        unsigned idx;
        unsigned cutdim;
        DataType cutval;
        meanSplit(rng, indices, count, idx, cutdim, cutval);

        node->DivDim = cutdim;
        node->DivVal = cutval;
        node->StartIdx = offset;
        node->EndIdx = offset + count;
        node->Lchild = divideTreeOnly(rng, indices, idx, offset);
        node->Rchild = divideTreeOnly(rng, indices+idx, count-idx, offset+idx);
      }

      return node;
    }
	 */

	void meanSplit(std::mt19937& rng, unsigned* indices, unsigned count, unsigned& index, unsigned& cutdim, DataType& cutval){
		size_t veclen_ = features_.get_cols();
		DataType* mean_ = new DataType[veclen_];
		DataType* var_ = new DataType[veclen_];
		memset(mean_,0,veclen_*sizeof(DataType));
		memset(var_,0,veclen_*sizeof(DataType));

		/* Compute mean values.  Only the first SAMPLE_NUM values need to be
          sampled to get a good estimate.
		 */
		unsigned cnt = std::min((unsigned)SAMPLE_NUM+1, count);
		for (unsigned j = 0; j < cnt; ++j) {
			const DataType* v = features_.get_row(indices[j]);
			for (size_t k=0; k<veclen_; ++k) {
				mean_[k] += v[k];
			}
		}
		DataType div_factor = DataType(1)/cnt;
		for (size_t k=0; k<veclen_; ++k) {
			mean_[k] *= div_factor;
		}

		/* Compute variances (no need to divide by count). */

		for (unsigned j = 0; j < cnt; ++j) {
			const DataType* v = features_.get_row(indices[j]);
			for (size_t k=0; k<veclen_; ++k) {
				DataType dist = v[k] - mean_[k];
				var_[k] += dist * dist;
			}
		}

		/* Select one of the highest variance indices at random. */
		cutdim = selectDivision(rng, var_);

		cutval = mean_[cutdim];

		unsigned lim1, lim2;

		planeSplit(indices, count, cutdim, cutval, lim1, lim2);
		//cut the subtree using the id which best balances the tree
		if (lim1>count/2) index = lim1;
		else if (lim2<count/2) index = lim2;
		else index = count/2;

		/* If either list is empty, it means that all remaining features
		 * are identical. Split in the middle to maintain a balanced tree.
		 */
		if ((lim1==count)||(lim2==0)) index = count/2;
		delete[] mean_;
		delete[] var_;
	}
	void planeSplit(unsigned* indices, unsigned count, unsigned cutdim, DataType cutval, unsigned& lim1, unsigned& lim2){
		/* Move vector indices for left subtree to front of list. */
		int left = 0;
		int right = count-1;
		for (;; ) {
			while (left<=right && features_.get_row(indices[left])[cutdim]<cutval) ++left;
			while (left<=right && features_.get_row(indices[right])[cutdim]>=cutval) --right;
			if (left>right) break;
			std::swap(indices[left], indices[right]); ++left; --right;
		}
		lim1 = left;//lim1 is the id of the leftmost point <= cutval
		right = count-1;
		for (;; ) {
			while (left<=right && features_.get_row(indices[left])[cutdim]<=cutval) ++left;
			while (left<=right && features_.get_row(indices[right])[cutdim]>cutval) --right;
			if (left>right) break;
			std::swap(indices[left], indices[right]); ++left; --right;
		}
		lim2 = left;//lim2 is the id of the leftmost point >cutval
	}
	int selectDivision(std::mt19937& rng, DataType* v){
		int num = 0;
		size_t topind[RAND_DIM];

		//Create a list of the indices of the top RAND_DIM values.
		for (size_t i = 0; i < features_.get_cols(); ++i) {
			if ((num < RAND_DIM)||(v[i] > v[topind[num-1]])) {
				// Put this element at end of topind.
				if (num < RAND_DIM) {
					topind[num++] = i;            // Add to list.
				}
				else {
					topind[num-1] = i;         // Replace last element.
				}
				// Bubble end value down to right location by repeated swapping. sort the varience in decrease order
				int j = num - 1;
				while (j > 0  &&  v[topind[j]] > v[topind[j-1]]) {
					std::swap(topind[j], topind[j-1]);
					--j;
				}
			}
		}
		// Select a random integer in range [0,num-1], and return that index.
		int rnd = rng()%num;
		return (int)topind[rnd];
	}
	void getMergeLevelNodeList(Node* node, size_t treeid, int deepth){
		if(node->Lchild != NULL && node->Rchild != NULL && deepth < ml){
			deepth++;
			getMergeLevelNodeList(node->Lchild, treeid, deepth);
			getMergeLevelNodeList(node->Rchild, treeid, deepth);
		}else if(deepth == ml){
			mlNodeList.push_back(std::make_pair(node,treeid));
		}else{
			error_flag = true;
			if(deepth < max_deepth)max_deepth = deepth;
		}
	}
	Node* SearchToLeaf(Node* node, size_t id){
		if(node->Lchild != NULL && node->Rchild !=NULL){
			if(features_.get_row(id)[node->DivDim] < node->DivVal)
				return SearchToLeaf(node->Lchild, id);
			else
				return SearchToLeaf(node->Rchild, id);
		}
		else
			return node;
	}int cc = 0;
	void mergeSubGraphs(size_t treeid, Node* node){
		if(node->Lchild != NULL && node->Rchild != NULL){
			mergeSubGraphs(treeid, node->Lchild);
			mergeSubGraphs(treeid, node->Rchild);

			size_t numL = node->Lchild->EndIdx - node->Lchild->StartIdx;
			size_t numR = node->Rchild->EndIdx - node->Rchild->StartIdx;
			size_t start,end;
			Node * root;
			if(numL < numR){
				root = node->Rchild;
				start = node->Lchild->StartIdx;
				end = node->Lchild->EndIdx;
			}else{
				root = node->Lchild;
				start = node->Rchild->StartIdx;
				end = node->Rchild->EndIdx;
			}

			for(;start < end; start++){

				size_t feature_id = LeafLists[treeid][start];

				Node* leaf = SearchToLeaf(root, feature_id);
				for(size_t i = leaf->StartIdx; i < leaf->EndIdx; i++){
					size_t tmpfea = LeafLists[treeid][i];
					DataType dist = distance_->compare(
							features_.get_row(tmpfea), features_.get_row(feature_id), features_.get_cols());

					{LockGuard g(*nhoods[tmpfea].lock);
					if(knn_graph[tmpfea].size() < params_.S || dist < knn_graph[tmpfea].begin()->distance){
						Candidate<DataType> c1(feature_id, dist);
						knn_graph[tmpfea].insert(c1);
						if(knn_graph[tmpfea].size() > params_.S)knn_graph[tmpfea].erase(knn_graph[tmpfea].begin());


					}
					else if(nhoods[tmpfea].nn_new.size() < params_.S * 2){

						nhoods[tmpfea].nn_new.push_back(feature_id);

					}
					}
					{LockGuard g(*nhoods[feature_id].lock);
					if(knn_graph[feature_id].size() < params_.S || dist < knn_graph[feature_id].begin()->distance){
						Candidate<DataType> c1(tmpfea, dist);
						knn_graph[feature_id].insert(c1);
						if(knn_graph[feature_id].size() > params_.S)knn_graph[feature_id].erase(knn_graph[feature_id].begin());

					}
					else if(nhoods[feature_id].nn_new.size() < params_.S * 2){

						nhoods[feature_id].nn_new.push_back(tmpfea);

					}
					}
				}
			}
		}
	}

	typedef std::set<Candidate<DataType>, std::greater<Candidate<DataType>> > CandidateHeap;


protected:
	enum
	{
		/**
		 * To improve efficiency, only SAMPLE_NUM random values are used to
		 * compute the mean and variance at each level when building a tree.
		 * A value of 100 seems to perform as well as using all values.
		 */
		SAMPLE_NUM = 100,
		/**
		 * Top random dimensions to consider
		 *
		 * When creating random trees, the dimension on which to subdivide is
		 * selected at random from among the top RAND_DIM dimensions with the
		 * highest variance.  A value of 5 works well.
		 */
		RAND_DIM=5
	};

	int TreeNum;
	int TreeNumBuild;
	int ml;   //merge_level
	int max_deepth;
	int veclen_;
	//DataType* var_;
	omp_lock_t rootlock;
	bool error_flag;
	//DataType* mean_;
	std::vector<Node*> tree_roots_;
	std::vector< std::pair<Node*,size_t> > mlNodeList;
	std::vector<std::vector<unsigned>> LeafLists;
	USING_BASECLASS_SYMBOLS

	//kgraph code

	static void GenRandom (std::mt19937& rng, unsigned *addr, unsigned size, unsigned N) {
		for (unsigned i = 0; i < size; ++i) {
			addr[i] = rng() % (N - size);
		}
		std::sort(addr, addr + size);
		for (unsigned i = 1; i < size; ++i) {
			if (addr[i] <= addr[i-1]) {
				addr[i] = addr[i-1] + 1;
			}
		}
		unsigned off = rng() % N;
		for (unsigned i = 0; i < size; ++i) {
			addr[i] = (addr[i] + off) % N;
		}
	}


	void DFSbuild(Node* node, std::mt19937& rng, unsigned* indices, unsigned count, unsigned offset){
		//omp_set_lock(&rootlock);
		//std::cout<<node->treeid<<":"<<offset<<":"<<count<<std::endl;
		//omp_unset_lock(&rootlock);
		if(count <= params_.TNS){
			node->DivDim = -1;
			node->Lchild = NULL;
			node->Rchild = NULL;
			node->StartIdx = offset;
			node->EndIdx = offset + count;
			//add points

		}else{
			unsigned idx;
			unsigned cutdim;
			DataType cutval;
			meanSplit(rng, indices, count, idx, cutdim, cutval);
			node->DivDim = cutdim;
			node->DivVal = cutval;
			node->StartIdx = offset;
			node->EndIdx = offset + count;
			Node* nodeL = new Node(); Node* nodeR = new Node();
			node->Lchild = nodeL;
			nodeL->treeid = node->treeid;
			DFSbuild(nodeL, rng, indices, idx, offset);
			node->Rchild = nodeR;
			nodeR->treeid = node->treeid;
			DFSbuild(nodeR, rng, indices+idx, count-idx, offset+idx);
		}
	}

	void DFStest(unsigned level, unsigned dim, Node* node){
		if(node->Lchild !=NULL){
			DFStest(++level, node->DivDim, node->Lchild);
			//if(level > 15)
			std::cout<<"dim: "<<node->DivDim<<"--cutval: "<<node->DivVal<<"--S: "<<node->StartIdx<<"--E: "<<node->EndIdx<<" TREE: "<<node->treeid<<std::endl;
			if(node->Lchild->Lchild ==NULL){
				std::vector<unsigned>& tmp = LeafLists[node->treeid];
				for(unsigned i = node->Rchild->StartIdx; i < node->Rchild->EndIdx; i++)
					std::cout<<features_.get_row(tmp[i])[node->DivDim]<<" ";
				std::cout<<std::endl;
			}
		}
		else if(node->Rchild !=NULL){
			DFStest(++level, node->DivDim, node->Rchild);
		}
		else{
			std::cout<<"dim: "<<dim<<std::endl;
			std::vector<unsigned>& tmp = LeafLists[node->treeid];
			for(unsigned i = node->StartIdx; i < node->EndIdx; i++)
				std::cout<<features_.get_row(tmp[i])[dim]<<" ";
			std::cout<<std::endl;
		}
	}
	void buildTrees(){
		unsigned N = features_.get_rows();
		unsigned seed = 1998;
		std::mt19937 rng(seed);
		nhoods.resize(N);
		g.resize(N);
		boost::dynamic_bitset<> visited(N, false);
		knn_graph.resize(N);
		for (auto &nhood: nhoods) {
			//nhood.nn_new.resize(params_.S * 2);
			nhood.pool.resize(params_.L+1);
			nhood.radius = std::numeric_limits<float>::max();
		}


		//build tree
		std::vector<int> indices(N);
		LeafLists.resize(TreeNum);
		std::vector<Node*> ActiveSet;
		std::vector<Node*> NewSet;
		for(unsigned i = 0; i < (unsigned)TreeNum; i++){
			Node* node = new Node;
			node->DivDim = -1;
			node->Lchild = NULL;
			node->Rchild = NULL;
			node->StartIdx = 0;
			node->EndIdx = N;
			node->treeid = i;
			tree_roots_.push_back(node);
			ActiveSet.push_back(node);
		}
#pragma omp parallel for
		for(unsigned i = 0; i < N; i++)indices[i] = i;
#pragma omp parallel for
		for(unsigned i = 0; i < (unsigned)TreeNum; i++){
			std::vector<unsigned>& myids = LeafLists[i];
			myids.resize(N);
			std::copy(indices.begin(), indices.end(),myids.begin());
			std::random_shuffle(myids.begin(), myids.end());
		}
		omp_init_lock(&rootlock);
		while(!ActiveSet.empty() && ActiveSet.size() < 1100){
#pragma omp parallel for
			for(unsigned i = 0; i < ActiveSet.size(); i++){
				Node* node = ActiveSet[i];
				unsigned mid;
				unsigned cutdim;
				DataType cutval;
				std::mt19937 rng(seed ^ omp_get_thread_num());
				std::vector<unsigned>& myids = LeafLists[node->treeid];

				meanSplit(rng, &myids[0]+node->StartIdx, node->EndIdx - node->StartIdx, mid, cutdim, cutval);

				node->DivDim = cutdim;
				node->DivVal = cutval;
				//node->StartIdx = offset;
				//node->EndIdx = offset + count;
				Node* nodeL = new Node(); Node* nodeR = new Node();
				nodeR->treeid = nodeL->treeid = node->treeid;
				nodeL->StartIdx = node->StartIdx;
				nodeL->EndIdx = node->StartIdx+mid;
				nodeR->StartIdx = nodeL->EndIdx;
				nodeR->EndIdx = node->EndIdx;
				node->Lchild = nodeL;
				node->Rchild = nodeR;
				omp_set_lock(&rootlock);
				if(mid>params_.S)NewSet.push_back(nodeL);
				if(nodeR->EndIdx - nodeR->StartIdx > params_.S)NewSet.push_back(nodeR);
				omp_unset_lock(&rootlock);
			}
			ActiveSet.resize(NewSet.size());
			std::copy(NewSet.begin(), NewSet.end(),ActiveSet.begin());
			NewSet.clear();
		}
#pragma omp parallel for
		for(unsigned i = 0; i < ActiveSet.size(); i++){
			Node* node = ActiveSet[i];
			//omp_set_lock(&rootlock);
			//std::cout<<i<<":"<<node->EndIdx-node->StartIdx<<std::endl;
			//omp_unset_lock(&rootlock);
			std::mt19937 rng(seed ^ omp_get_thread_num());
			std::vector<unsigned>& myids = LeafLists[node->treeid];
			DFSbuild(node, rng, &myids[0]+node->StartIdx, node->EndIdx-node->StartIdx, node->StartIdx);
		}
	}
    void outputVisitBucketNum(){}

	void initGraph(){
		//initial
		unsigned N = features_.get_rows();
		unsigned seed = 1998;
		std::mt19937 rng(seed);
		nhoods.resize(N);
		g.resize(N);
		boost::dynamic_bitset<> visited(N, false);
		knn_graph.resize(N);
		for (auto &nhood: nhoods) {
			//nhood.nn_new.resize(params_.S * 2);
			nhood.pool.resize(params_.L+1);
			nhood.radius = std::numeric_limits<float>::max();
		}


		//build tree
		std::vector<int> indices(N);
		LeafLists.resize(TreeNum);
		std::vector<Node*> ActiveSet;
		std::vector<Node*> NewSet;
		for(unsigned i = 0; i < (unsigned)TreeNum; i++){
			Node* node = new Node;
			node->DivDim = -1;
			node->Lchild = NULL;
			node->Rchild = NULL;
			node->StartIdx = 0;
			node->EndIdx = N;
			node->treeid = i;
			tree_roots_.push_back(node);
			ActiveSet.push_back(node);
		}
#pragma omp parallel for
		for(unsigned i = 0; i < N; i++)indices[i] = i;
#pragma omp parallel for
		for(unsigned i = 0; i < (unsigned)TreeNum; i++){
			std::vector<unsigned>& myids = LeafLists[i];
			myids.resize(N);
			std::copy(indices.begin(), indices.end(),myids.begin());
			std::random_shuffle(myids.begin(), myids.end());
		}
		omp_init_lock(&rootlock);
		while(!ActiveSet.empty() && ActiveSet.size() < 1100){
#pragma omp parallel for
			for(unsigned i = 0; i < ActiveSet.size(); i++){
				Node* node = ActiveSet[i];
				unsigned mid;
				unsigned cutdim;
				DataType cutval;
				std::mt19937 rng(seed ^ omp_get_thread_num());
				std::vector<unsigned>& myids = LeafLists[node->treeid];

				meanSplit(rng, &myids[0]+node->StartIdx, node->EndIdx - node->StartIdx, mid, cutdim, cutval);

				node->DivDim = cutdim;
				node->DivVal = cutval;
				//node->StartIdx = offset;
				//node->EndIdx = offset + count;
				Node* nodeL = new Node(); Node* nodeR = new Node();
				nodeR->treeid = nodeL->treeid = node->treeid;
				nodeL->StartIdx = node->StartIdx;
				nodeL->EndIdx = node->StartIdx+mid;
				nodeR->StartIdx = nodeL->EndIdx;
				nodeR->EndIdx = node->EndIdx;
				node->Lchild = nodeL;
				node->Rchild = nodeR;
				omp_set_lock(&rootlock);
				if(mid>params_.S)NewSet.push_back(nodeL);
				if(nodeR->EndIdx - nodeR->StartIdx > params_.S)NewSet.push_back(nodeR);
				omp_unset_lock(&rootlock);
			}
			ActiveSet.resize(NewSet.size());
			std::copy(NewSet.begin(), NewSet.end(),ActiveSet.begin());
			NewSet.clear();
		}
#pragma omp parallel for
		for(unsigned i = 0; i < ActiveSet.size(); i++){
			Node* node = ActiveSet[i];
			//omp_set_lock(&rootlock);
			//std::cout<<i<<":"<<node->EndIdx-node->StartIdx<<std::endl;
			//omp_unset_lock(&rootlock);
			std::mt19937 rng(seed ^ omp_get_thread_num());
			std::vector<unsigned>& myids = LeafLists[node->treeid];
			DFSbuild(node, rng, &myids[0]+node->StartIdx, node->EndIdx-node->StartIdx, node->StartIdx);
		}
		//DFStest(0,0,tree_roots_[0]);
		//build tree completed

		for(size_t i = 0; i < (unsigned)TreeNumBuild; i++){
			getMergeLevelNodeList(tree_roots_[i], i ,0);
		}

#pragma omp parallel for	
		for(size_t i = 0; i < mlNodeList.size(); i++){
			mergeSubGraphs(mlNodeList[i].second, mlNodeList[i].first);
		}


#pragma omp parallel
		{
#ifdef _OPENMP
			std::mt19937 rng(seed ^ omp_get_thread_num());
#else
			std::mt19937 rng(seed);
#endif
			std::vector<unsigned> random(params_.S + 1);

#pragma omp for
			for (unsigned n = 0; n < N; ++n) {
				auto &nhood = nhoods[n];
				Points &pool = nhood.pool;
				if(nhood.nn_new.size()<params_.S*2){
					nhood.nn_new.resize(params_.S*2);
					GenRandom(rng, &nhood.nn_new[0], nhood.nn_new.size(), N);
				}


				GenRandom(rng, &random[0], random.size(), N);
				nhood.L = params_.S;
				nhood.Range = params_.S;
				while(knn_graph[n].size() < params_.S){
					unsigned rand_id = rng() % N;
					DataType dist = distance_->compare(
							features_.get_row(n), features_.get_row(rand_id), features_.get_cols());
					Candidate<DataType> c(rand_id,dist);
					knn_graph[n].insert(c);
				}

				//omp_set_lock(&rootlock);
				//if(knn_graph[n].size() < nhood.L)std::cout<<n<<":"<<knn_graph[n].size()<<std::endl;
				//omp_unset_lock(&rootlock);
				unsigned i = 0;
				typename CandidateHeap::reverse_iterator it = knn_graph[n].rbegin();
				for (unsigned l = 0; l < nhood.L; ++l) {
					if (random[i] == n) ++i;
					auto &nn = nhood.pool[l];
					nn.id = it->row_id;//random[i++];
					nhood.nn_new[l] = it->row_id;
					nn.dist = it->distance;//distance_->compare(features_.get_row(n), features_.get_row(nn.id), features_.get_cols());
					nn.flag = true;it++;
					//if(it == knn_graph[n].rend())break;
				}
				sort(pool.begin(), pool.begin() + nhood.L);
			}
		}
		knn_graph.clear();
#ifdef INFO
		std::cout<<"initial completed"<<std::endl;
#endif
	}

};

}
#endif

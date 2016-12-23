// Copyright (C) 2016 Cong Fu <731097343@qq.com>. All Rights Reserved.
#ifndef EFANNA
#define EFANNA
#include "general/distance.hpp"
#include "general/matrix.hpp"
#include "general/params.hpp"
#include "algorithm/init_indices.hpp"
namespace efanna{
template <typename DataType>
class FIndex{
public:
	typedef InitIndex<DataType> IndexType;

	FIndex(const Matrix<DataType>& features, Distance<DataType>* d, const IndexParams& params)
	: index_params_(params)
	{
		init_algorithm init_index_type= params.init_index_type;
		index_params_ = params;
		initIndex_ = create_index_by_type(init_index_type, features, params, d);
	}

	virtual ~FIndex () {
	}
	void buildIndex()
	{
		initIndex_->buildIndex();
	}
	void buildTrees()
	{
		initIndex_->buildTrees();
	}
	void knnSearch(int k, const Matrix<DataType>& query){
		initIndex_->knnSearch(k, query);
	}
	void saveIndex(char* filename){
		initIndex_->saveIndex(filename);
	}
	void loadIndex(char* filename){
		initIndex_->loadIndex(filename);
	}
	void saveTrees(char* filename){
		initIndex_->saveTrees(filename);
	}
	void loadTrees(char* filename){
		initIndex_->loadTrees(filename);
	}
	void saveGraph(char* filename){
		initIndex_->saveGraph(filename);
	}
	void loadGraph(char* filename){
		initIndex_->loadGraph(filename);
	}
	void saveResults(char* filename){
		initIndex_->saveResults(filename);
	}
	void setSearchParams(int epochs, int init_num, int extend_to, int search_trees = 0, int search_lv=-1, int search_method = 0){
		initIndex_->setSearchParams(epochs, init_num, extend_to,search_trees, search_lv, search_method);
	}
	size_t getGraphSize(){
		return initIndex_->getGraphSize();
	}
	std::vector<unsigned> getGraphRow(unsigned row_id){
		return initIndex_->getGraphRow(row_id);
	}
	void outputVisitBucketNum(){
		initIndex_->outputVisitBucketNum();
	}

private:
	/** Pointer to actual index class */
	IndexType* initIndex_;
	/** Parameters passed to the index */
	IndexParams index_params_;

};
}
#endif

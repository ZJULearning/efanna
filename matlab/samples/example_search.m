% add the path of efanna library
addpath('../')

% use fvecs_read util function provided to load datasets and queries
dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
query = fvecs_read('~/data/sift/sift_query.fvecs');
disp('Query size:');
disp(size(query));
% params: dataset, number of trees for overall building, conquer-to-depth, number of iterations, L, checkK, K, S, number of trees for graph building
% all params should be strictly corresponds to the params in building process of loaded graph and trees
% click following link for more information: https://github.com/fc731097343/efanna/blob/master/README.md
ef = efanna(dataset, 16, 8, 8, 30, 25, 10, 10, 8);
% build, or load graph & trees
ef.load_trees('./sift.trees');
ef.load_graph('./sift.graph');
%%% ef.build_index();

% search params: number of trees to use, number of epoches, pool size factor, extend factor, searching methods
% more details can be found in README.md of efanna root as well.
ef.set_search_params(16, 4, 1200, 200, 0);
% params: required number of returned neighbors (i.e. k for knn, here are searching for 10-nn)
ef.knn_search(10, query);
ef.save_result('./sift_search.result');

% add the path of efanna library
addpath('../')

% use fvecs_read util function provided to load datasets and queries
dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
query = fvecs_read('~/data/sift/sift_query.fvecs');
disp('Query size:');
disp(size(query));
% params: dataset, index name, distance name, rnn_used, trees, mlevel, epoces, L, checkK, K, S, build_trees
% click following link for more information: https://github.com/fc731097343/efanna/blob/master/README.md
ef = efanna(dataset, 'kdtreeub', 'l2', true, 16, 8, 8, 30, 25, 10, 10);
% build, or load graph & trees
%% ef.load_trees('./sift.trees');
%% ef.load_graph('./sift.graph');
ef.build_index();

% search params: search_trees, search_depth, search_epoch, search_extend, pool_size
% more details can be found in README.md of efanna root as well.
ef.set_search_params(16, 11, 4, 100, 100);
% params: k for knn (here are searching for 1-nn)
ef.knn_search(1, query);
ef.save_result('./sift_search.result');

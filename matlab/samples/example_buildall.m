% add the path of efanna library
addpath('../')

% use fvecs_read util function provided to load datasets
dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
% params: dataset, number of trees for overall building, conquer-to-depth, number of iterations, L, checkK, K, S, number of trees for graph building
% click following link for more information: https://github.com/fc731097343/efanna/blob/master/README.md
ef = efanna(dataset, 16, 8, 8, 30, 25, 10, 10, 8);
% build graph and get a sparse matrix describing the NN results
spmat = ef.build_index();
disp('Adjacency matrix of KNN graph acquired. Shape:');
disp(size(spmat));
disp('Number of non-zero elements in adjacency matrix of KNN graph:');
disp(nnz(spmat));
% save trees and graph
ef.save_trees('./sift.trees');
ef.save_graph('./sift.graph');

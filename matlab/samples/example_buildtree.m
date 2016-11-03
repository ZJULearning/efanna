% add the path of efanna library
addpath('../')

% use fvecs_read util function provided to load datasets
dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
% params: dataset, index name, distance name, rnn_used, trees, mlevel, epoces, L, checkK, K, S
% click following link for more information: https://github.com/fc731097343/efanna/blob/master/README.md
ef = efanna(dataset, 'kdtreeub', 'l2', true, 16);
% build trees and save
ef.build_trees();
ef.save_trees('./sift.trees');

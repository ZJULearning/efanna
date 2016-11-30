% add the path of efanna library
addpath('../')

% use fvecs_read util function provided to load datasets
dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
% params: dataset, number of trees to build
% click following link for more information: https://github.com/fc731097343/efanna/blob/master/README.md
ef = efanna(dataset, 16);
% build trees and save
ef.build_trees();
ef.save_trees('./sift.trees');

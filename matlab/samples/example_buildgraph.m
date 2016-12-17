% add the path of efanna library
addpath('../')

% use fvecs_read util function provided to load datasets
dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
% params: dataset, number of trees for graph building, conquer-to-depth, number of iterations, L, checkK, K, S
% click following link for more information: https://github.com/fc731097343/efanna/blob/master/README.md
ef = efanna(dataset, 8, 8, 7, 30, 25, 10, 10);
% build graph and get a sparse matrix describing the NN results
% the sparse matrix is row-wise organized, i.e. the first row lies the 0-1 vector describes k nearest neighbours for node 1.
tic;
spmat = ef.build_index();
toc;

gdgraph = ivecs_read('~/data/sift/sift_10NN_groundtruth.graph');
spmat = spmat';
[row,col]=find(spmat);
row = row-1;
nCorrect = 0;
for i=1:1000000
    for j=1:10
        if(find(gdgraph(:,i)==row((i-1)*10+j)))
            nCorrect = nCorrect + 1;
        end
    end
end 
disp(['10NN accuracy: ',num2str(nCorrect/10000000)]);
% save graph
ef.save_graph('./sift.graph');

addpath('../')

dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
ef = efanna(dataset, 'kdtreeub', 'l2', true, 16, 5, 8, 25, 30, 10);
spmat = ef.build_index();
disp('Adjacency matrix of KNN graph acquired. Shape:');
disp(size(spmat));
disp('Non-zero elements in adjacency matrix of KNN graph:');
disp(nnz(spmat));

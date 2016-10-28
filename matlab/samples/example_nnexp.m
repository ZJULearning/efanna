addpath('../')

dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
ef = efanna(dataset, 'nnexp', 'l2', 10, 10);
ef.build_index();
ef.save_index('~/data/sift/sift_nnexp.idx');


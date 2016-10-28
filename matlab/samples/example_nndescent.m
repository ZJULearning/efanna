addpath('../')

dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
ef = efanna(dataset, 'nndescent', 'l2', true, 10, 25, 30);
ef.build_index();
ef.save_index('~/data/sift/sift_nndescent.idx');


addpath('../')

dataset = fvecs_read('~/data/sift/sift_base.fvecs');
disp('Data size:');
disp(size(dataset));
ef = efanna(dataset, 'kdtreeub', 'l2', true, 16, 5, 8, 25, 30, 10);
ef.build_index();
ef.save_graph('~/data/sift/sift_kdtree.graph');

query = fvecs_read('~/data/sift/sift_query.fvecs');
ef.set_search_params(500, 5, 11, 16, 500);
ef.knn_search(1, query);
ef.save_result('~/data/sift/sift_kbtreeub.result');


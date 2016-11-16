classdef efanna < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle;
        index_name;
    end
    methods
        % Constructor
        function this = efanna(data, index_name, dist_name, varargin)
            this.objectHandle = findex('new', data, index_name, dist_name, varargin{:});
            this.index_name = index_name;
        end
        % Destructor
        function delete(this)
            findex('destruct', this.objectHandle);
        end
        
        function graph_mat = build_index(this)
            graph_mat = findex('build_index', this.objectHandle);
        end

        function build_trees(this)
            findex('build_trees', this.objectHandle);
        end

        function load_index(this, filename)
            findex('load_index', this.objectHandle, filename);
        end

        function save_index(this, filename)
            findex('save_index', this.objectHandle, filename);
        end

        function graph_mat = load_graph(this, filename)
            graph_mat = findex('load_graph', this.objectHandle, filename);
        end

        function save_graph(this, filename)
            findex('save_graph', this.objectHandle, filename);
        end

        function load_trees(this, filename)
            findex('load_trees', this.objectHandle, filename);
        end

        function save_trees(this, filename)
            findex('save_trees', this.objectHandle, filename);
        end
        % set search params (different index may share the same knn_search with different params)
        function set_search_params(this, varargin)
            findex('set_search_params', this.objectHandle, this.index_name, varargin{:});
        end

        function knn_search(this, k, query)
            findex('knn_search', this.objectHandle, k, query);
        end

        function save_result(this, filename)
            findex('save_result', this.objectHandle, filename);
        end
    end
end

classdef Dijkstra
    
    properties
        numnodes            % number of detected nodes 
        type = 'min'        % if using min (distances) or max (probabilities)
        start_nodes         % start nodes
        end_nodes           % end nodes 
        weights             % weight graph 
        nodemap             % indices that enables extraction of values of a matrix at the original node indices 
        revnodemap          % vector of node indices to help revert back to original to get peaks  
        
        dist_graph          % distance graph (unweighted)
        path_len            % how long paths should be
        graph
        hweight  = 2        % weight going horizontally

        thresh_matrix 
        
        %%% more for debugging if anything 
        rowidx      % esp helpful when we have the shortest path in term of node numbers and we want to see 
                            % which row in the original matrix each node
                            % belongs to (b/c we know that each step in the
                            % trajectory should correspond to a unique
                            % cluster 
        colidx      % if done correctly using the path as an index to colidx should yield 1:nclusters
        exttype             % extrema type peaks (1) or valley (0)
    end
    
    methods
        

        function obj = Dijkstra(varargin)
        % Initialization of the Dijkstra class, setting any properties if
        % specified
            if ~isempty(varargin)
                for arg = 1:length(varargin)
                    if isprop(obj, varargin{arg})
                        obj = setval(obj, varargin{arg}, varargin{arg+1});
                    end                    
                end                
            end 
        end
        
        
        function  obj = setval(obj, str, input)
        % set the value of the property if applicable 
        
            if isprop(obj, str)
                %checked_input = check(obj, str, input);
                checked_input = input;
                obj.(str) = checked_input;
            else
                error([str, ' is not a property of Dijkstra class'])
            end
            
        end
        
        function obj = pred_nodes(obj, nodepred_matrix, thresh)
        % find nodes that we consider are probable enough to be a peak
        % candidate
        % Inputs:
        %       obj         [Dijkstra class] 
        %       nodepred_matrix  [NxM]  probability matrix with probabilities 
        %                                        of a peak Ni being in a
        %                                        certain cluster Mi
        %       thresh:          [dbl]  threshold to reduce nodepred_matrix

            % eliminate nodes that have an extremely low threshold 
            nodepred_matrix(nodepred_matrix < thresh) = 0;
            
            % set the thresh holded matrix for easier viewing later 
            obj = obj.setval('thresh_matrix', nodepred_matrix);
            
            % confirm there is at least 1 non-zero components in each gaussian
            nnzero_cand = sum(nodepred_matrix ~= 0, 1);
            if sum(nnzero_cand ~= 0) ~= length(nnzero_cand)
                disp('Number of non-zero elements in each cluster')
                disp(nnzero_cand)
                error('some clusters have no candidates: ')

            end

            % calculate the number of starting and ending candidates
            nstart_cand = sum(nodepred_matrix(:, 1) ~= 0);
            nend_cand = sum(nodepred_matrix(:, end)~= 0);
            
            % get the linear indices of non-zero values (nodes)
            revnodemap_temp = find(nodepred_matrix ~=0);
            
            % save the reverse mapping from nodes back to a matrix 
            obj = obj.setval('revnodemap', revnodemap_temp);
            % save the number of nodes detected 
            obj = obj.setval('numnodes', length(obj.revnodemap));
            
            % create a node list 
            nodelist = [1:obj.numnodes]';
            
            % calculate the original row and column of each node to the
            % original matrix so we can convert back
            [rowidx_temp, colidx_temp] = ind2sub(size(nodepred_matrix), obj.revnodemap);
            
            % make row and column indices into column vectors 
            if size(rowidx_temp, 1) < size(rowidx_temp, 2)
               rowidx_temp = rowidx_temp'; 
            end
            if size(colidx_temp, 1) < size(colidx_temp, 2)
               colidx_temp = colidx_temp'; 
            end
            
            % we save these because they are easier to parse
            obj = obj.setval('rowidx', rowidx_temp);
            obj = obj.setval('colidx', colidx_temp);
            
            % now convert starting candidates and ending candidates to their node index
            obj = obj.setval('start_nodes', nodelist(1:nstart_cand));
            obj = obj.setval('end_nodes', nodelist(end-nend_cand + 1:end));

            % ideally don't have to do more than 20 nodes (but this can be altered) ....
            if length(obj.end_nodes) * length(obj.start_nodes) > 20
                error('too many candidate pairings required, consider reducing threshold')
            end            
        end
        function obj = gen_weight_multi(obj)
        % generate weight vector to encourage/discourage certain node paths
        % ok now we have a set of rules here we asssume columns are clusters and
        % rows are points. We also assume that clusters that refer to peaks
        % and clusters that refer to valleys always alternate, Similarly
        % rows that refer to extrema also alternate peak, valley, peak
        % Ex: Node assignment example 
        %       | c1 (peak)    c2 (valley)    c3 (peak)     c4 (valley)
        %       |______________________________________________________
        %       |
        %    p1 |   1               7            13             19              
        %    v1 |   2               8            14             20
        %    p2 |   3               9            15             21
        %    v2 |   4              10            16             22
        %    p3 |   5              11            17             23
        %    v3 |   6              12            18             24
        % 
        % Distance Weight Assignment (source row, source col) -> (dest row, dest col) 
        %   
        %   1. col_i -> col_i:  0 
        %               (peaks/valleys from the same gaussian are never connected)
        %   2. sign(row_k) = a -> , sign(row_k+n) = ~a & col_i -> col_i+1 ; n > 0 : 1
        %               (we cant go from a peak in a peak cluster to a peak in a valley cluster)
        %   3. col_i, row_j -> col_i+1, row_j-n; n > 0: 4*value
        %               (this assumes that the next cluster chooses a 
        %               peak BEFORE the peak chosen for the current cluster, we also 
        %               might want to consider just doing n = 1 and setting n > 1 
        %               to 0)
        %   4. col_i -> col_i-1 : 0
        %               (can't travel backwards in clusters)
        %   5. col_i -> col_i+n; n>1 : 0
        %               (we cannot skip clusters)
        %   5. col_i, row_i -> col_i+1, row_i+1: value
        % Inputs:
        %       obj         [Dijkstra class] 
        
            % set the source row (ascending delays), columns (ascending
            % clusters), 
            source_row = repmat(obj.rowidx, 1, obj.numnodes);
            source_col= repmat(obj.colidx, 1, obj.numnodes);
            source_sign = obj.exttype(source_row);
            
            dest_row = source_row';
            dest_col = source_col';
            dest_sign = obj.exttype(dest_row);
            
            % determine what distance metric to use (d2, gnpost_prob, ...)
            % for distances we want to penalize large distances for weights
            % 
            if strcmp(obj.type, 'min')
                exp_fac = 1;
            % for probabilities 
            else 
                exp_fac = -1;
            end
            
            % set up the placeholder for weights (holds emphasis for distances)
            weights_temp = zeros(size(source_row));


            % rule 2 (peak -> valley or valley -> peak) & (col_i ->
            % col_i+1) & (row_h - row_k > 0) (lower diagonal)
            rsame_rule = (dest_sign - source_sign ~= 0) & (dest_col - source_col == 1) ...
                & (dest_row - source_row > 0);
            weights_temp(rsame_rule) = 1; % look at this the regular results were 4^exp_fac

            % rule 3 (peak -> valley or valley -> peak) & (col_i ->
            % col_i+1) & (row_h - row_k = -1) (upper diagonal)
            rlwer_rule = (dest_sign - source_sign ~= 0) & (dest_col - source_col == 1) ...
                & (dest_row - source_row == -1);
            weights_temp(rlwer_rule) = obj.hweight^(exp_fac);

            % rule 5 (uber upper diagonal)
            cuppr_rule = (dest_sign - source_sign ~= 0) & (dest_col - source_col == 1) ...
                & (dest_row - source_row < -1);
            weights_temp(cuppr_rule) = obj.hweight^(2*exp_fac);

            
            obj = obj.setval('weights', weights_temp);
        end
        
        function obj = gen_weight(obj)
        % generate weight vector to encourage/discourage certain node paths
        % ok now we have a set of rules here we asssume columns are clusters and
        % rows are points 
        %   1. col_i -> col_i:  0 
        %               (peaks from the same gaussian are never connected)
        %   2. col_i, row_j -> col_i+1, row_j: 2*value 
        %               (we discourage the same peak being chosen across clusters)
        %   3. col_i, row_j -> col_i+1, row_j-n; n > 0: 4*value
        %               (this assumes that the next cluster chooses a 
        %               peak BEFORE the peak chosen for the current cluster, we also 
        %               might want to consider just doing n = 1 and setting n > 1 
        %               to 0)
        %   4. col_i -> col_i-1 : 0
        %               (can't travel backwards in clusters)
        %   5. col_i -> col_i+n; n>1 : 0
        %               (we cannot skip clusters)
        %   5. col_i, row_i -> col_i+1, row_i+1: value
        % Inputs:
        %       obj         [Dijkstra class] 
            source_row = repmat(obj.rowidx, 1, obj.numnodes);
            source_col= repmat(obj.colidx, 1, obj.numnodes);
            if sum(obj.thresh_matrix) == 5
                obj;
            end
            dest_row = source_row';
            dest_col = source_col';
            
            % determine what distance metric to use (d2, gnpost_prob, ...)
            % for distances we want to penalize large distances for weights
            % 
            if strcmp(obj.type, 'min')
                exp_fac = 1;
            % for probabilities 
            else 
                exp_fac = -1;
            end
            % set up the placeholder for weights (holds emphasis for distances)
            weights_temp = zeros(size(source_row));


            % rules 1 and 4 
            col_rule = (dest_col - source_col <= 0);
            weights_temp(col_rule) = 0;

            % rule 2
            rsame_rule = (dest_row - source_row == 0 & dest_col - source_col == 1);
            weights_temp(rsame_rule) = obj.hweight^exp_fac; % look at this the regular results were 4^exp_fac

            % rule 3 
            rlwer_rule = (dest_row - source_row < 0 &dest_col - source_col > 0);
            weights_temp(rlwer_rule) = 2*obj.hweight^(exp_fac);

            % rule 5
            cuppr_rule = (dest_col - source_col > 1);
            weights_temp(cuppr_rule) = 0;

            % rule 6
            pos_rule  = (dest_col - source_col == 1 & dest_row - source_row > 0);
            weights_temp(pos_rule) = 1;
            
            obj = obj.setval('weights', weights_temp);
        end
        
        
        function obj = prep_dijkstra(obj, nodepred_matrix, thresh, dist_matrix)
            % here we prep all the things necessary for running the
            % dijkstra algorithm. This means 
            % 1. mapping of gmm output -> trajectory 
            % 2. weight matrix 
            % We take in threshold to threshold matrix and nodepred_matrix
            % which allows us to predict what hte nodes are 
            
            obj = obj.setval('path_len', size(nodepred_matrix, 2));
            
            %%%%%%% 
            % -- Find potential nodes and get the starting and ending nodes
            %%%%%%%
            obj = obj.pred_nodes(nodepred_matrix, thresh);
            
            % create a node list 
            nodelist = [1:obj.numnodes]';
            
            % just create a matrix where source contains the first node and dest
            % contains the 2nd node the source will travel to. Do the same for the
            % column and rows associated with the nodes. Now we can see that 
            % source(2, 1) -> dest(2, 1) = 2 -> 1; in other words from node 2 to node 1
            source = repmat(nodelist, 1, obj.numnodes);
            dest = source';

            
            % save the indices to convert any matrix to node form 
            nodemap_temp = obj.revnodemap(dest);

            % 
            obj = obj.setval('nodemap', nodemap_temp);
            
            %%%%%%% 
            % -- Calculate the weighting matrix 
            %%%%%%%
            obj = obj.gen_weight_multi();
            
            %%%%%%% 
            % -- Generate initial trajectory matrix (unweighted)
            %%%%%%%            
            obj = obj.setval('dist_graph', dist_matrix(obj.nodemap));
            
            % calculate the 
            graph = obj.dist_graph .* obj.weights; 
            obj = obj.setval('graph', graph);
        end
        function W = dijkstra(obj, graph, source)
        %---------------------------------------------------
        % Dijkstra Algorithm
        % author : Dimas Aryo
        % email : mr.dimasaryo@gmail.com
        % This code has been modified slightly to suit the purposes of fid
        % point detection in SCG
        % Inputs
        %       obj         [Dijkstra class] 
        %       graph       [NxN]  Distance from every permutation of node
        %                               to node
        %       source      [dbl]  Node to start
        %---------------------------------------------------       
        graph(graph == 0) = inf;
        graph= obj.exchangenode(graph,1,source);

        lengthA=size(graph,1);
        W=zeros(lengthA);
        for i=2 : lengthA
            W(1,i)=i;
            W(2,i)=graph(1,i);
        end

        for i=1 : lengthA
            D(i,1)=graph(1,i);
            D(i,2)=i;
        end

        D2=D(2:length(D),:);
        L=2;
        while L<=(size(W,1)-1)
            L=L+1;
            D2=sortrows(D2,1);
            k=D2(1,2);
            W(L,1)=k;
            D2(1,:)=[];
            for i=1 : size(D2,1)
                if D(D2(i,2),1)>(D(k,1)+graph(k,D2(i,2)))
                    D(D2(i,2),1) = D(k,1)+graph(k,D2(i,2));
                    D2(i,1) = D(D2(i,2),1);
                end
            end

            for i=2 : length(graph)
                W(L,i)=D(i,1);
            end
        end
        end
        function [L, e] = listdijkstra(obj, L,W,s,d)
        %---------------------------------------------------
        % Dijkstra Algorithm
        % author : Dimas Aryo
        % email : mr.dimasaryo@gmail.com
        % Find the path in term of node numbers from one node to another
        % This code has been modified slightly to suit the purposes of fid
        % point detection in SCG
        % Inputs
        %       obj         [Dijkstra class] 
        %       L      [NxN]  Distance from every permutation of node
        %                               to node
        %       source      [dbl]  Node to start
        %--------------------------------------------------- 
            e=W(size(W,1),d);
            index=size(W,1);
            while index>0
                if W(2,d)==W(size(W,1),d)
                    L=[L s];
                    index=0;
                else
                    index2=size(W,1);
                    while index2>0
                        if W(index2,d)<W(index2-1,d)
                            if W(index2,1)==s
                                L = [L 1];
                            else
                                L=[L W(index2,1)];
                            end
                            L =obj.listdijkstra(L,W,s,W(index2,1));
                            index2=0;
                        else
                            index2=index2-1;
                        end
                        index=0;
                    end
                end
            end
        end
        function [route_dist, node_routes, node_dist_routes] = calc_paths(obj, start_nodes, end_nodes)
            

            
            route_dist = zeros(length(start_nodes)*length(end_nodes), 1);
            node_routes = zeros(length(start_nodes)*length(end_nodes), obj.path_len);
            node_dist_routes = zeros(length(start_nodes)*length(end_nodes), obj.path_len);
            for i = 1:length(start_nodes)
                W = obj.dijkstra(obj.graph, start_nodes(i));
                for j = 1:length(end_nodes)
                    [flipped_route, route_dist(i*j)] = obj.listdijkstra([end_nodes(j)],W,start_nodes(i),end_nodes(j));
%                     if length(flipped_route) ~= obj.path_len
%                         node_routes(i*j, :) = [];
%                         node_dist_routes(i*j, :) = [];
%                         route_dist(i*j) = [];
%                     else
                        node_routes(i*j, :) = flip(flipped_route);
                        starts = node_routes(i*j, 1:end-1);
                        ends = node_routes(i*j, 2:end);
                        node_dist_routes(i*j, :) = [0, diag(obj.graph(starts, ends))'];
%                     end
                end
            end
            % note sure why some distances are zero? 
            bad_routes = find(route_dist == 0);
            node_routes(bad_routes, :) = [];
            route_dist(bad_routes) = [];
        end

        
    end
    
    methods (Static)
        
        function G = exchangenode(G,a,b)
        % enables calculating distances if b is not actually the first node
        % (a)
        
            %Exchange element at column a with element at column b;
            buffer=G(:,a);
            G(:,a)=G(:,b);
            G(:,b)=buffer;

            %Exchange element at row a with element at row b;
            buffer=G(a,:);
            G(a,:)=G(b,:);
            G(b,:)=buffer;
            
        end
        

    end
end
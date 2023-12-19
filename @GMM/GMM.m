 classdef GMM 
    
    properties 
        nclusters 
        ndim
        mu          % should be a nclusters x ndim matrix
        Sigma           % should be a ndim x ndim x ncluster matrix 
        ComponentProportion           % should be a 1 x ncluster matrix 
        maxiter = 100                 % number of iterations to cap
        model 
        SharedCovariance = 0
        CovarianceType = 'diagonal'
    end    
        
    methods 
        
        function obj = GMM(varargin)
            if ~isempty(varargin)
                for arg = 1:length(varargin)
                    if isprop(obj, varargin{arg})
                        obj = setval(obj, varargin{arg}, varargin{arg+1});
                    end                    
                end                
            end 
        end
        
        function obj = setStart(obj, data, groups)
        % Takes the data, parses it using the groups label, and gets the 
        % mean and covariance for each label. We also assume the priors 
        % for each label is uniform 
        % Input:
        %       data    [Nxndim]    dataset to warm start for GMM 
        %       groups  [Nx1]       labels for each data point 
        % 

            ndim = size(data, 2);
                        
            if ndim ~= obj.ndim
               error('data should be a Nxndim matrix') 
            end
            
            % initalize covariances, means, and priors 
            mu = zeros(obj.nclusters, ndim);
            if strcmp(obj.CovarianceType, 'diagonal')
                Sigma = repmat(eye(ndim), 1, 1, obj.nclusters);
            else
                Sigma = repmat(eye(ndim), 1, 1, obj.nclusters);
            end
            ComponentProportion = ones(1, obj.nclusters)*1/obj.nclusters;

            % for each label find the covariances and means 
            for cluster_i = 1:obj.nclusters
                mu(cluster_i, :) = mean(data(groups == cluster_i, :), 1); 
                Sigma_i = cov(data(groups == cluster_i, :));
                if strcmp(obj.CovarianceType, 'diagonal')
                    Sigma(:, :, cluster_i) = diag((diag(Sigma_i)));
                else
                    Sigma(:, :, cluster_i) = Sigma_i;
                end
            end

            % set the means, covariances, and priors for the GMM 
            obj = setval(obj, 'mu', mu);
            obj = setval(obj, 'Sigma', Sigma);
            obj = setval(obj, 'ComponentProportion', ComponentProportion);
            
        end            
        
        function  obj = setval(obj, str, input)
            if isprop(obj, str)
                
                checked_input = check(obj, str, input);
                obj.(str) = checked_input;
            else
                error([str, ' is not a property of GMM class'])
            end
        end
                  
        function obj = fitModel(obj, data)
            S = struct;
            S.mu = obj.mu;
            if strcmp(obj.CovarianceType, 'diagonal')
                S.Sigma = zeros(1, obj.ndim, obj.nclusters);
                for i = 1:obj.nclusters
                    S.Sigma(:, :, i) = diag(obj.Sigma(:, :, i))';                 
                end
            end

            
            
            S.ComponentProportion = obj.ComponentProportion;
            obj.model = fitgmdist(data, obj.nclusters, 'Start', S, ...
                'SharedCovariance', logical(obj.SharedCovariance), 'CovarianceType', obj.CovarianceType);
            
            if strcmp(obj.CovarianceType, 'full')
                obj.mu = obj.model.mu;
                for i = 1:obj.nclusters
                    obj.Sigma(:, :, i) = diag(obj.model.Sigma(:, :, i));
                end
            end
            
        end
        
        function [idx, post_pdf, npost_pdf, mdist] = evalPoints(obj, data, varargin)
            % nll not necessary
            % probs is stupid af it just takes the logpdf subtracts the max
            % value - basically making the most likely prob 0 - and then reverses
            % the log by taking the exp... 
            %[idx, nll, probs, logpdf, d2] = cluster(obj.model, data);
            %post_pdf = zeros(size(probs));
            
            % type 
            %       'gprob':    uses cond prob P(cluster_i | x_i)
            %       'mdist':    uses malhabonis dist
            type = 'gprob'; 
            delay = 0;
            if ~isempty(varargin)
                for arg = 1:length(varargin)
                    if strcmp('type', varargin{arg}); type = varargin{arg+1};end
                    if strcmp('delay', varargin{arg}); delay = varargin{arg+1}; end
                end
            end
               
            % multi dim and want to calculate vlaues using all dimensions
            if obj.ndim ~= 1 && delay == 0
%                 if strcmp(obj.model.CovarianceType, 'full')
%                     CovType = 2;
%                 else
%                     CovType = 1;
%                 end
%                 [log_lh, mdist] = obj.wdensity_edit(data, obj.model.mu, obj.model.Sigma, ...
%                     obj.model.ComponentProportion, obj.model.SharedCovariance, CovType);
%                 [ll, post, logpdf]= obj.estep(log_lh);
%                 post_pdf = exp(log_lh);
%                 npost_pdf = post_pdf./sum(post_pdf, 2);
                post_pdf = zeros(size(data, 1), 1);
                mdist = zeros(size(data, 1), 1);
                
                for i = 1:obj.nclusters
                    post_pdf(:, i) = mvnpdf(data, obj.mu(i, :), obj.Sigma(:, :, i))*obj.model.ComponentProportion(i);
                    mdist(:, i) = diag(sqrt((data - obj.mu(i, :))*inv(diag(diag(obj.Sigma(:, :, i))))*(data - obj.mu(i, :))'));
                end
                % conditional prob p(g|x)
                npost_pdf = post_pdf./sum(post_pdf, 2);
                
            % multi dim but just want to use one dim 
            elseif obj.ndim ~= 1 && delay > 0
                
                % conditional prob p(x|g) and mahalabonis distance
                post_pdf = zeros(size(data, 1), 1);
                mdist = zeros(size(data, 1), 1);

                % get the gmm vars and means associated with delay dim
                Sigma = squeeze(obj.Sigma(delay, delay, :));
                mu = obj.mu(:, delay);

                % for each cluster
                for i = 1:obj.nclusters
                    % calculate conditional prob and mahalbonis dist
                    post_pdf(:, i) = normpdf(data(:, delay), mu(i), Sigma(i))*obj.model.ComponentProportion(i);
                    mdist(:, i) = sqrt((data(:, delay) - mu(i)).^2/Sigma(i));
                end
                
                % conditional prob p(g|x)
                npost_pdf = post_pdf./sum(post_pdf, 2);
                
            
            % if we are using 1 dimension or if we have 2 and just want ot
            % evaluate on delay parameters 
            else
                % conditional prob p(x|g)
                post_pdf = zeros(size(data, 1), 1);
                mdist = zeros(size(data, 1), 1);
                if obj.ndim == 1
                    Sigma = obj.Sigma;
                    mu = obj.mu;
                else
                    Sigma = squeeze(obj.Sigma(1, 1, :));
                    mu = obj.mu(:, 1);
                end
                for i = 1:obj.nclusters
                    post_pdf(:, i) = normpdf(data(:, 1), mu(i), Sigma(i))*obj.model.ComponentProportion(i);
                    mdist(:, i) = sqrt((data(:, 1) - mu(i)).^2/Sigma(i));
                end
                
                % conditional prob p(g|x)
                npost_pdf = post_pdf./sum(post_pdf, 2);
                
            end
            
            % finds the most probable cluster for each sample using either
            % the post posterior or the malahabonis distance 
            if strcmp(type, 'gprob')
                idx = obj.findCluster(post_pdf, 'max');
            else
                idx = obj.findCluster(mdist, 'min');
            end

        end
        

        function genGMMcontour(obj, data, varargin)
        % visualization function 
            circle = 0;
            type = 'positive';
            if ~isempty(varargin)
                for arg = 1:length(varargin)
                   if strcmp('points', varargin{arg}); points = varargin{arg+1}; circle = 1; end  
                   if strcmp('type', varargin{arg}); type = varargin{arg+1};; end  
                end
            end
            % only do this if ndim = 2 or 3
            [N, ndim1] = size(data);
            if obj.ndim ~= ndim1
                error('data should have same dimension as gmm')
            end
            if obj.ndim > 3 
                error('plotting capabilities only valid to 2 or 3 dimensions')
            end
            %figure
            pointsz = 20;
            if obj.ndim == 1
                
                gmPDF = @(x) arrayfun(@(x0) pdf(obj.model,[x0 ]), x);
                scatter(data, gmPDF(data)); hold on;
                plotdim = [min(data), max(data)];
                fplot(gmPDF,plotdim)
                
                if circle 
                    scatter(points, gmPDF(points), [], 'r');
                end
            elseif obj.ndim == 2

                gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(obj.model,[x0 y0]),x,y);
                scatter(data(:,1),data(:,2),pointsz,'.') 
                hold on
                % https://www.mathworks.com/help/stats/gmdistribution.pdf.html
                plotdim = [min(data(:, 1)), max(data(:, 1)), min(data(:, 2)), max(data(:, 2))];
                fcontour(gmPDF,plotdim)
                xlabel('Dim 1')
                ylabel('Dim 2')
            elseif obj.ndim == 3
                 scatter3(data(:,1),data(:,2), data(:, 3),pointsz,'.') 
%                 hold on
%                 gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(obj.model,[x0 y0]),x,y);
%                 plotdim = [min(data(:, 1)), max(data(:, 1)), min(data(:, 2)), max(data(:, 2))];
%                 fcontour(gmPDF,plotdim)
            end
                
        end      
        
        function output = check(obj, str, input)
                output = input;
                
                % means should be a ndim x ngauss
                if strcmp(str, 'means')

                    if ~ischar(input)
                        [nclusters1, ndim1] = size(input);
                        if ndim1 ~= obj.ndim | nclusters1 ~= obj.nclusters
                            error(['means is %d x %d but should be a %d x %d matrix (nclusters x ndim)'], ...
                                nclusters1, ndim1, obj.nclusters, obj.dim)
                        end
                    elseif strcmp(input, 'default')
                        output = zeros(obj.ndim, obj.nclusters);
                    else 
                        error('input should be a matrix or "default"')
                    end
                
                    
                % sigma should be a ndim x ndim x ncluster
                % matrices
                elseif strcmp(str, 'sigma')
                    if ~ischar(input)
                        [ndim1, ndim2, nclusters1] = size(input);
                        if ndim1 ~= obj.ndim | ndim2 ~= obj.ndim | nclusters1 ~= obj.nclusters
                            error(['sigma is %d x %d x %d but should be a %d x %d x %d matrix (ndim x ndim x nclusters)'], ...
                                ndim1, ndim2, nclusters1, obj.ndim, obj.ndim, obj.nclusters)
                        end

                        
                    elseif strcmp(input, 'default')
                        
                        output = zeros(obj.ndim, obj.ndim, obj.nclusters);                        
                        for i = 1:obj.nclusters
                            output(:, :, i) = eye(obj.ndim);
                        end
                        
                    else 
                        error('input should be a cell vector or "default"')
                    end         
                
                
                % prior should be a 1x nclusters vector
                elseif strcmp(str, 'prior')

                    if ~ischar(input)
                        [ndim1, nclusters1] = size(input);
                        if ndim1 ~= 1 | nclusters1 ~= obj.nclusters
                            error(['prior is %d x %d but should be a 1 x %d matrix (ndim x ndim x nclusters)'], ...
                                ndim1, nclusters1, obj.ndim, obj.nclusters)
                        end
                        if sum(input) ~= 1
                            error('input should sum up to 1')
                        end
                    elseif strcmp(input, 'default')
                        output = 1/obj.nclusters * ones(obj.nclusters, 1);
                    else 
                        error('input should be a vector or "default"')
                    end
                    
                end
        end
    end
    
    methods (Static)
        function  [ll, post, logpdf]=estep(log_lh,prob_th)
        %ESTEP E-STEP for Gaussian mixture distribution
        %   LL = ESTEP(LOG_LH) returns the loglikelihood of data in LL.  LOG_LH
        %   is the log of component conditional density weighted by the component
        %   probability.
        %
        %   [LL, POST] = ESTEP(LOG_LH) returns the posterior probability in the
        %   matrix POST. POST(i,j) is the posterior  probability of point i
        %   belonging to cluster j. 
        %
        %   [LL, POST] = ESTEP(LOG_LH,PROB_TH) set the probablities that is smaller
        %   than PROB_TH to zero.
        %
        %   [LL, POST, DENSITY] = ESTEP(LOG_LH) returns the logs of the pdf values
        %   of data in the vector density.
        %
        %   Copyright 2007-2016 The MathWorks, Inc.

        maxll = max(log_lh,[],2);
        %minus maxll to avoid underflow
        post = exp(log_lh-maxll);
        %density(i) is \sum_j \alpha_j P(x_i| \theta_j)/ exp(maxll(i))
        density = sum(post,2);
        logpdf = log(density) + maxll;
        ll = sum(logpdf); 
        post = post./density;%normalize posteriors

        %Set small posteriors to zero for efficiency
        %Currently, the following steps are only performed in the fitting phase
        if nargin > 1
            post(post<(prob_th))=0;
            density = sum(post,2);
            post = post./density;%renormalize posteriors
        end
        end
        
        function [log_lh,mahalaD] = wdensity_edit(X, mu, Sigma, p, sharedCov, CovType)
               log_prior = log(p);
                [n,d]=size(X);
                k=size(mu,1);
                log_lh = zeros(n,k,'like',X);
                if nargout > 1
                  mahalaD = zeros(n,k,'like',X);
                end
                logDetSigma = -Inf;
                for j = 1:k
                    if sharedCov
                        if j == 1
                            if CovType == 2 % full covariance
                                [L,f] = chol(Sigma);
                                diagL = diag(L);
                                if (f ~= 0)|| any(abs(diagL) < eps(max(abs(diagL)))*size(L,1))
                                    error(message('stats:gmdistribution:wdensity:IllCondCov'));
                                end
                                logDetSigma = 2*sum(log(diagL));
                            else %diagonal
                                L = sqrt(Sigma);
                                if  any(L < eps( max(L))*d)
                                      error(message('stats:gmdistribution:wdensity:IllCondCov'));
                                end
                                logDetSigma = sum( log(Sigma) );
                            end
                        end
                    else %different covariance
                        if CovType == 2 %full covariacne
                            % compute the log determinant of covariance
                            [L,f] = chol(Sigma(:,:,j) );
                            diagL = diag(L);
                            if (f ~= 0) || any(abs(diagL) < eps(max(abs(diagL)))*size(L,1))
                                 error(message('stats:gmdistribution:wdensity:IllCondCov'));
                            end
                            logDetSigma = 2*sum(log(diagL));
                        else %diagonal covariance
                            L = sqrt(Sigma(:,:,j)); % a vector
                            if  any(L < eps(max(L))*d)
                                 error(message('stats:gmdistribution:wdensity:IllCondCov'));
                            end
                            logDetSigma = sum(log(Sigma(:,:,j)) );
                        end
                    end

                    if CovType == 2
                         log_lh(:,j) = sum(((X - mu(j,:))/L).^2, 2); 
                    else %diagonal covariance
                         log_lh(:,j) = sum(((X - mu(j,:))./L).^2, 2); 
                    end

                    if nargout > 1
                         mahalaD(:,j) = log_lh(:,j);
                    end
                    log_lh(:,j) = -0.5*(log_lh(:,j) + logDetSigma);
                end
                %log_lh is a N by K matrix, log_lh(i,j) is log \alpha_j(x_i|\theta_j)
                log_lh = log_lh + log_prior - d*log(2*pi)/2;
        end  
        
        function idx = findCluster(value, type)
            % assumes value is a N x ncluster matrix and finds for each row
            % the best cluster based on a simple min/max search along each
            % row 
            % type determines if it's going to be min or max 
            if strcmp(type, 'max')
                [~, linidx] = max(value, [], 2, 'linear');
                [~, idx] = ind2sub(size(value), linidx);
            else
                [~, linidx] = min(value, [], 2, 'linear');
                [~, idx] = ind2sub(size(value), linidx);  
            end
        end
        
        function [idx] = eval_labels(est_labels, true_labels)
        end
        
        
    end
    
end
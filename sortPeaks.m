function [sortedPeaks, varargout] = sortPeaks(scg_beats, fs, varargin)
% Function
% Let default be the  
% Input (required)
%       scg_beats    [MxN]      N beats of M samples each 
%       fs           [dbl]      sampling frequency 
%
% Inputs (optional)
%       Update flg   [flg]      update GMM parameter flag 
%       Dijktra flg  [flg]      Use Dijkstra for peak sorting flag 
%       NumClusters  [dbl]      number of clusters to track
%       Ndim         [dbl]      number of dimensions/features for each
%                                   cluster (1: peak delay, 2: peak amplitude, 3:
%                                   peak width)
%       Dijktresh    [dbl]      Threshold for dijkstra to reduce node candidates
%       Ncand        [dbl]      Max number of peaks to search for in each beat
%       Brange       [1x2 dbl]  starting and ending idx of beats to use in training (if applicable)
%       Axrange      [1x2 dbl]  starting and ending time (ms) to search for
%                                   relevant aortic peaks
%       Multikalman  [flg]      use multikalman implementation 
%       Singlekalman [flg]      use single kalman implementation 
%       Verbose      [flg]      Output text flag 
%       Plot         [flg]      Plot relevant plots flag 
%       Qconst       [dbl]      constant to scale covariance matrix for
%                                   kalman filter 
%       Influence    [dbl]      how much each peak will influence the gmm
%                                   update
%       HWeight      [dbl]      how much to penalize paths where the same
%                                   peak is selected for multiple clusters 
%       Extrema      [flg]      Track peaks and valleys flag 
%       Training     [KxM]          assume the argument provided is a clean K x M dataset
%       Manual       [Bx2]     first column is index, second column is
%                                   whether it is a peak or valley (1/0) 
% Notes 
%       This function will probably NOT work if scg_beats or the training
%       set provided contains NaNs
% Usage: 
% use multi kalman filter with dijkstra (default) 
% locs = sortPeaks(sepbeats_cln, fs, 'NClusters', 4, 'AXRange', [1, 250]);
% use single kalman filter w/o dijkstra or update
% locs = sortPeaks(sepbeats_cln, fs, 'NClusters', 4, 'AXRange', [1, 250], 'SingleKalman', 1, ...
%    'MultiKalman', 0, 'Dijkstra', 1, 'Update', 0);
% use multi kalman but penalize the same node being chosen for different clusters
% locs = sortPeaks(sepbeats_cln, fs, 'NClusters', 4, 'AXRange', [1, 250], 'Hweight', 4);
% use multi kalman but reduce influence of peaks on gmm update
% locs = sortPeaks(sepbeats_cln, fs, 'NClusters', 4, 'AXRange', [1, 250], 'Influence', 1/1000);
%srtpeakstoc = toc

% set the default arguments 
options = struct('model', 0, ...        % default model provided
    'dijkstra', 1, ...                  % default using dijkstra  
    'nclusters', 4, ...                 % default num clusters 
    'ndim', 1, ... 
    'dijkthresh', 0.1, ...              % threshold to consider for nodes 
    'training', 0, ...                  % training data provided 
    'ncand', 10, ...                    % number of candidate peaks 
    'brange', [1, 50], ...
    'axrange', [1, 250], ...            % aortic regions range (ms)
    'nst', 30, ...                      % default template size
    'err', 30, ...                      % default error 
    'multikalman', 1, ...               % algorithm defaults to peak;mu state vector
    'singlekalman', 0, ...              % algorithm can be set to peak state vector
    'verbose', 1, ...                   % set f0r debugging 
    'plot', 1, ...                      % set for debugging 
    'qconst', 0.2, ...                  % process noise covariance scalar
    'rconst', 1, ...                    % observation noise covariance scalar
    'influence' , [], ...               % how much the new peak will influence the overall mean 
    'hweight', 2, ...
    'custom', [], ...
    'ensavg', 0, ...
    'prealter', 0, ...
    'promconst', 10, ...                % the prominence constant 
    'del', 0, ...                       % 0: when using ndim > 1 calculate gmm with 2 dim 1: use gmm with time dim
    'relweights', 0, ...                % relweights for cluster relative informaiton 
    'mkfadaptive', 0, ...               % use dijk weighted multi-weighted adaptive 
    'dijkweight', 0, ...                % use dijk altered malhabonis weights 
    'strict', 0, ...                    % remove any crossings in training data
    'pkadj', 1, ...                     %
    'findextrema', 0, ...
    'mixd2', 0, ...                 % use a combination of distance and prob for R matrix
    'dijk_mixd2', 0, ...            % use a combination of distance and prob for dijk matrix
    'qmatrix_abs', 0, ...           % use the actual peak values rather than the differences for q matrix construction
    'gmmupdate', 1);                    % gmm update 

% get the field names 
optionNames = fieldnames(options);


% parse Name Value Pairs if provided
if ~isempty(varargin)
    for arg = 1:length(varargin)
        % check if a Name argument is given 
        if isstr(varargin{arg})
            
            % make case insenstive 
            field = lower(varargin{arg});
            
            % if optional name argument is valid
            if any(strcmp(field, optionNames))
                
                % for name value pairs with non-boolean values save tho
                if strcmp(field, 'gmm'); pos_gmm = varargin{arg+1};
                % otherwise save into the options struct; 
                elseif strcmp(field, 'training'); train_scg = varargin{arg+1}; options.training = 1;
                elseif strcmp(field, 'ensavg'); M = varargin{arg+1}; options.ensavg = 1; 
                else options.(field) = varargin{arg+1}; end
              
            % otherwise throw an error 
            else
                error([field, ' is not a valid Name Pair'])
            end
        end
    end
end

% if speicifed, apply exponential moving average on the beats
if options.ensavg
    scg_beats = cardio.general.ema(scg_beats, M, false);
end

% confirm that argument combination 
if options.multikalman && ~options.gmmupdate
    error('If Multikalman is used, GMMUpdate must be set to 1')
end

if options.singlekalman && options.multikalman 
    error('If SingleKalman specified, MultiKalman must be set to 0')
end

if options.mkfadaptive && (~options.dijkstra || ~options.multikalman)
    error('If MultiKalman Dijkstra Weighted Update is specified, Dijkstra must be set to 1, and MultiKalman must be set to 1')
end

% if options.pkadj && (~options.multikalman)
%     error('If peak adjustment multikalman must be set to 1')    
% end

% if speicfiied indicate the arguments used 
if options.verbose
    
    fprintf('Running sortPeaks with ')
    
    if options.dijkstra; fprintf('Dijkstra Sorting, ')
    elseif ~options.dijkstra; fprintf('No Dijkstra Sorting, '); end
        
    if options.multikalman; fprintf('MultiKalman');
        if options.mkfadaptive; fprintf(' with Dijkstra Weighted Updates, ')
        else fprintf(', '); end
    elseif options.singlekalman; fprintf('SingleKalman, ');            
    else; fprintf('No Kalman, '); end
    
    if options.relweights; fprintf('Relative GMM Distance Weights, ')
    elseif ~options.relweights; fprintf('No Relative GMM Distance Weights, '); end
    
    if ~options.prealter; fprintf('Kalmn filter multi using u_t = u_t + a*u_t-1, ')
    elseif options.prealter; fprintf('Kalman filter multi using u_t = u_t-1, '); end
    
    if options.multikalman; fprintf('Kalman Tracked GMM Update');
    elseif options.gmmupdate; fprintf('Single Step GMM Update');
    else fprintf('No GMM Update'); end
    
    fprintf('\n')
        
end   

% get the number of beats
nbeats = size(scg_beats, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       
%                   1a    Prepping Training Data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if no custom training set is provided 
if ~options.training
    
    if options.verbose
        fprintf('Using beats %d to %d for training\n', options.brange(1), options.brange(2));
    end

    % get the range of beats from the given dataset to extract
    ranges = options.brange(1):options.brange(2);
    
    % extract the training beats from testing dataset
    train_scg = scg_beats(:, ranges);

% if a custom training set is provided 
else
    
    if options.verbose
        fprintf('Using Training Data provided fro training\n')
    end
    
end

% determine the number of beats provided 
ntrain_beats = size(train_scg, 2);

% calculate influence value (how much the means of each cluster will change)
if isempty(options.influence)
    % if unsepecified use the number of training beats
    mu_val = 1/ntrain_beats;
else
    % otherwise use provided value
    mu_val = options.influence;
end

%if manual cluster locations were provided
if ~isempty(options.custom)
    
    % find approximate peak and valley features of aortic complexes using
    % sliding template
    [train_maxlocs, train_minlocs] = general.slideTemplate(train_scg, ...
        options.nst, options.err, fs, 'Closest', 'Axrange', options.axrange, ...
        'NumFeatures', options.nclusters, 'Manual', options.custom);
    
    % calculate total number of clusters
    options.nclusters = size(train_maxlocs, 1) + size(train_minlocs, 1);

    
% if no manual cluster locations were provided
else
    
    
    if options.findextrema
        template = mean(train_scg, 2);
        perc_lim = 0.3;
        [locs_seq, loc_signs] = find_extremas(template, options.axrange, fs, perc_lim, options.plot);
        options.custom = [locs_seq, loc_signs];
        [train_maxlocs, train_minlocs] = general.slideTemplate(train_scg, ...
            options.nst, options.err, fs, 'Closest', 'Axrange', options.axrange, ...
            'Manual', options.custom);
        % calculate total number of clusters
        options.nclusters = size(train_maxlocs, 1) + size(train_minlocs, 1);
       
    else
    % find approximate peak and valley features of aortic complexes using
    % sliding template with warm start of relevant features 
    [train_maxlocs, train_minlocs] = general.slideTemplate(train_scg, ...
        options.nst, options.err, fs, 'Closest', 'Axrange', options.axrange, ...
        'NumFeatures', options.nclusters, 'Emd');
    end
end

npkclusters = size(train_maxlocs, 1);
nvalclusters = size(train_minlocs, 1);

% save the peak locations and sort the training peaks 
train_maxlocs = train_maxlocs';
[~, sortidx] = sort(train_maxlocs(1, :), 'ascend');
train_maxlocs = train_maxlocs(:, sortidx);

% save the valley locations and sort the valley peaks 
train_minlocs = train_minlocs';
[~, sortidx] = sort(train_minlocs(1, :), 'ascend');
train_minlocs = train_minlocs(:, sortidx);

% concatenate peak/valley mean locations and their associated signs (1/0)
extremas = [mean(train_maxlocs), mean(train_minlocs)];
extrema_signs = [ones(1, size(train_maxlocs, 2)), zeros(1, size(train_minlocs, 2))];

% sort the peak and valley clusters and apply sorting on their signs and
% index number 
[~, extrema_gmm_ord] = sort(extremas, 'ascend');
extrema_signs = extrema_signs(extrema_gmm_ord);

% assign index number to peaks/valleys and apply sorting 
extrema_cols = [1:size(train_maxlocs, 2), 1:size(train_minlocs, 2)];
extrema_cols = extrema_cols(extrema_gmm_ord);

% if the ordering of signs is not peak, valley, peak or valley, peak,
% valley
if sum(diff(extrema_signs) == 0) ~= 0
    
    % find the locations where the peak/valley sequence fails
    poor_values = find(diff(extrema_signs) == 0);
    poor_values = poor_values + 1;
    
    % determine if the error is due to a peak or valley cluster
    pos_remove = extrema_cols(poor_values(extrema_signs(poor_values) == 1));
    neg_remove = extrema_cols(poor_values(extrema_signs(poor_values) == 0));
    
    % remove the peak or valley cluster that is causing the error and
    % update the number of peak/valley clusters
    if ~isempty(pos_remove)
        train_maxlocs(:, pos_remove) = [];
        npkclusters = size(train_maxlocs, 2);
    end
    if ~isempty(neg_remove)
        train_minlocs(:, neg_remove) = [];
        nvalclusters = size(train_minlocs, 2);
    end   
    
    % inidcate to the user that the number of clusters specified was lowered
    disp(['Was instructed to find ', num2str(options.nclusters), '. However only ', ...
        num2str(nvalclusters+npkclusters), ' were tracked'])
    
    % update total number of clusters
    options.nclusters = npkclusters + nvalclusters;
    
    % concatentate peak/valley mean locations and signs (1/0)
    extremas = [mean(train_maxlocs), mean(train_minlocs)];
    extrema_signs = [ones(1, size(train_maxlocs, 2)), zeros(1, size(train_minlocs, 2))];
        
    % sort peak/valley clusters and apply sorting on signs 
    [~, extrema_gmm_ord] = sort(extremas, 'ascend');
    extrema_signs = extrema_signs(extrema_gmm_ord);
    
    % assign index number to peaks/valleys and apply sorting 
    extrema_cols = [1:size(train_maxlocs, 2), 1:size(train_minlocs, 2)];
    extrema_cols = extrema_cols(extrema_gmm_ord);    
    
end

% if strictness of training set is specified (no crossings between clusters)
if options.strict
    
    % concatenate peaks/valley and sort based on sign order
    all_extremas = [train_maxlocs, train_minlocs];
    all_extremas = all_extremas(:, extrema_gmm_ord);
    
    % find beat indices where the locations are not strictly ascending
    crossings = sum((diff(all_extremas, [], 2) < 0) > 0, 2);
    
    % remove those beats and separate them back into peak and valley clusters
    all_extremas = all_extremas(~crossings, :);
    train_maxlocs = all_extremas(:, extrema_signs == 1);
    train_minlocs = all_extremas(:, extrema_signs == 0); 
    
    % recreate training set and update number of beats
    train_scg = train_scg(:, ~crossings);
    ntrain_beats = size(train_scg, 2);
    ranges = 1:ntrain_beats;
    
end

% plot the training data with peak and valley clusters identified
if options.plot
    
    figure;
    
    % plot the training data with peak and valley clusters identified
    subplot(2, 1, 1); 
    hold on; 
    
    % get the peak and valley amplitudes given peak/valley locations
    train_amps = general.indexMatrix(train_scg, train_maxlocs);
    train_minamps = general.indexMatrix(train_scg, train_minlocs);
    
    % concatenate peaks/valley locations and sort based on sign order
    all_extremas = [train_maxlocs, train_minlocs];
    all_extremas = all_extremas(:, extrema_gmm_ord);
    
    % concatenate peaks/valley amplitudes and sort based on sign order
    all_amps = [train_amps, train_minamps];
    all_amps = all_amps(:, extrema_gmm_ord);
    
    % plot the locations and amplitude paris 
    for i = 1:options.nclusters
        scatter(all_extremas(:, i), all_amps(:, i));
    end
    
    % plot the training signal set
    plot(train_scg); 
    title('Training SCGs with Fid Point estimated')
    
    % for color consistency get the plots associated with the scatter plot
    % and set them to the first children to use default colors
    h = get(gca, 'Children');
    set(gca, 'Children', [h(end-options.nclusters+1:end); h(1:end-options.nclusters)])
    
    % plot just the locations 
    subplot(2, 1, 2);
    plot(all_extremas); xlabel('Beats')
    title('Fid Points')
    ylabel('delay (samples)')
    
end

% if for some reason (mainly because of custom inputs) the signs do not
% flip indicate they should to the user
if sum(diff(extrema_signs) == 0) ~= 0
   error('Tracked peaks should be switching from peak and valleys') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       
%                1b       Training GMM
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a gmm object for the peak clusters
pos_gmm = GMM('nclusters', npkclusters, 'ndim', options.ndim);

% create a dataset format for the GMM to use 
train_posdata = zeros([size(train_maxlocs), options.ndim]);

% first dimension is always peak location 
train_posdata(:, :, 1) = train_maxlocs;

% if more than 1 dimension is specified for training
if options.ndim > 1
    
    % get the prominences and widths from dataset
    [pos_cand, neg_cand, pos_prom, neg_prom] = ...
        findAO(train_scg, 20, fs, [options.axrange(1), options.axrange(2)], 1);
    
    % set placeholder for prominences of the training peak locations
    train_maxproms = zeros(size(train_maxlocs));
    
    % for each peak cluster 
    for i = 1:npkclusters
        
        % get the indices of the peaks that best match the training peak locations
        prom_idx = abs(train_maxlocs(:, i) - pos_cand) == min(abs(train_maxlocs(:, i)-pos_cand), [], 2);
        
        % find the beats where multiple candidates were found and find
        % first instance candidate
        dup_rows = find(sum(prom_idx, 2) > 1);
        [~, good_idx] = max(prom_idx(dup_rows, :), [], 2);
        
        % reorganize such that each beat only has 1 candidate 
        prom_idx(dup_rows, :) = zeros(length(dup_rows), size(prom_idx, 2));
        prom_idx(sub2ind(size(prom_idx), dup_rows, good_idx)) = 1;
        
        % index and take all relevant prominence values 
        train_maxproms(:, i) = pos_prom(prom_idx);
    end
    
    % save the prominences in the 2nd dimension 
    train_posdata(:, :, 2) = train_maxproms;
    
end

% save training data and prep training data and labels to be used by gmm 
unprep_train_posdata = train_posdata;
[train_posdata, train_poslabels] = prep_data(train_posdata, 1:ntrain_beats);

% warm start training by finding initial means/sigmas from training data
pos_gmm = pos_gmm.setStart(train_posdata, train_poslabels);
pos_gmm = pos_gmm.fitModel(train_posdata);


% create a gmm object for vallye clusters
neg_gmm = GMM('nclusters', nvalclusters, 'ndim', options.ndim);

train_negdata = zeros([size(train_minlocs), options.ndim]);

% set the data to the first dimensions 
train_negdata(:, :, 1) = train_minlocs;

% for each vallye cluster
if options.ndim > 1
    
    % set placeholder for prominences of the training valley locations
    train_minproms = zeros(size(train_minlocs));
    
    % for each valley cluster
    for i = 1:size(train_minlocs, 2)
        
        % get the indices of the valleys that best match the training valley locations
        prom_idx = abs(train_minlocs(:, i) - neg_cand) == min(abs(train_minlocs(:, i)-neg_cand), [], 2);
        
        % find the beats where multiple candidates were found and find
        % first instance candidate
        dup_rows = find(sum(prom_idx, 2) > 1);
        [~, good_idx] = max(prom_idx(dup_rows, :), [], 2);
        
        % reorganize such that each beat only has 1 candidate 
        prom_idx(dup_rows, :) = zeros(length(dup_rows), size(prom_idx, 2));
        prom_idx(sub2ind(size(prom_idx), dup_rows, good_idx)) = 1;
        
        % index and take all relevant prominence values
        train_minproms(:, i) = neg_prom(prom_idx);
    end 

    % save the prominences in the 2nd dimension 
    train_negdata(:, :, 2) = train_minproms;
    
end

% save training data and prep training data and labels to be used by gmm 
unprep_train_negdata = train_negdata;
[train_negdata, train_neglabels] = prep_data(train_negdata, 1:ntrain_beats);

% warm start training by finding initial means/sigmas from training data
neg_gmm = neg_gmm.setStart(train_negdata, train_neglabels);
neg_gmm = neg_gmm.fitModel(train_negdata);



% if using multi dimensions 
if options.ndim ~= 1
    %neg_gmm.Sigma(end, end, :) = neg_gmm.Sigma(end, end, :)*options.promconst;
    %pos_gmm.Sigma(end, end, :) = pos_gmm.Sigma(end, end, :)*options.promconst;
    
    neg_gmm.Sigma = neg_gmm.Sigma*options.promconst;
    pos_gmm.Sigma = pos_gmm.Sigma*options.promconst;
end

% if sepcified show plot of distributions found with training data labeled
if options.plot
    figure;
    neg_gmm.genGMMcontour(train_negdata)
    pos_gmm.genGMMcontour(train_posdata)
end

% find the order of the peak/valley clusters and their associated signs (0/1)
[~, extrema_gmm_ord] = sort([pos_gmm.mu(:, 1); neg_gmm.mu(:, 1)], 'ascend');
extrema_signs = [ones(size(pos_gmm.mu(:, 1))); zeros(size(neg_gmm.mu(:, 1)))];
extrema_signs = extrema_signs(extrema_gmm_ord);

% determine if the first extrema is a vally or peak
extrema_first = extrema_signs(1);

% extract whether calculating distance just on time or all available
% dimensions 
del = options.del;

% initialize placeholder for relative distances of peaks to other clusters
unprep_pos_relations = zeros(options.nclusters, npkclusters);

% for each peak cluster
for i = 1:npkclusters
    
    % calculate distances of peaks using just the time axis with peak cluster i to both peak and valley gmms
    [~, ~, ~, pos_d2] = pos_gmm.evalPoints(permute(unprep_train_posdata(:, i, :), [1, 3, 2]), 'type', 'mdist', 1);
    [~, ~, ~, neg_d2] = neg_gmm.evalPoints(permute(unprep_train_posdata(:, i, :), [1, 3, 2]), 'type', 'mdist', 1);
    
    % merge distances and sort based on gmm cluster location order
    d2 = [pos_d2, neg_d2];
    d2 = d2(:, extrema_gmm_ord);
    
    % calculate average distances between all gmms and all valleys in valley cluster i
    dist_pos_relations = mean(d2, 1)';
    
    %%%%%%%%%%%%%%%% I WILL NEED THIS %%%%%%%%%%%%%%%%%%%%%
%     [~, idx] = min(dist_pos_relations);
%     sign_relations = ones(size(dist_pos_relations));
%     sign_relations(1:idx-1) = -1;
%     unprep_pos_relations(:, i) = sign_relations .* dist_pos_relations;

    unprep_pos_relations(:, i) = dist_pos_relations;
end

% initialize placeholder for realtive distances with other clusters
unprep_neg_relations = zeros(options.nclusters, nvalclusters);
for i = 1:size(unprep_train_negdata, 2)
    
    % calculate distances using just time time axis of valleys with valley cluster i to both peak and valley gmms
    [~, ~, ~, pos_d2] = pos_gmm.evalPoints(permute(unprep_train_negdata(:, i, :), [1, 3, 2]), 'type', 'mdist', 'delay', 1);
    [~, ~, ~, neg_d2] = neg_gmm.evalPoints(permute(unprep_train_negdata(:, i, :), [1, 3, 2]), 'type', 'mdist', 'delay', 1);
    
    % merge distances and sort based on gmm order 
    d2 = [pos_d2, neg_d2];
    d2 = d2(:, extrema_gmm_ord);
    
    % calculate average distances between all gmms and all valleys in valley cluster i
    dist_neg_relations = mean(d2, 1)';
    
    %%%%%%%%%%%%%%%% I will need this for 
    %%%%%%%%%%%%%%%%%% Consider using the x/y dimensions of d2 
%     [~, idx] = min(dist_neg_relations);
%     sign_relations = ones(size(dist_neg_relations));
%     sign_relations(1:idx-1) = -1;
%    unprep_neg_relations(:, i) = sign_relations .* dist_neg_relations;
    unprep_neg_relations(:, i) = dist_neg_relations;
end

% this will hold the tracked peak indices for the data provided 
sortedPeaks = zeros(nbeats, options.nclusters);
sortedPeaks_d2 = zeros(nbeats, options.nclusters);
sortedPeaks_probs = zeros(nbeats, options.nclusters);
sortedPeaks_nprobs = zeros(nbeats, options.nclusters);
allmu = sortedPeaks;

% info for each datapoint will have the location, width, d2, gnprob
info = struct([]);

sortedProms = sortedPeaks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          1c. Initialize Kalman Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.multikalman
    % here we assume state_vector = [p1, p2, p3, p4, u1, u2, u3, u4]
    nstates = options.nclusters*2; % based on number of kalman filters prior ]  
   
    % predicted big boi
    P_pred = zeros(nstates, nstates);

    % residual big boi 
    res = zeros(nbeats, nstates);

    % state estimation matrix 
    stateEst = zeros(nbeats, nstates);   

    % create state transition matrix (assume mu_t = mu_t-1) 
    ss_A = eye(nstates);
    
    % if specified our state vector is now mu_t+1 = mu_t-1 + a*(mu_t-1 + peak_t-1) 
    if ~options.prealter
        ss_A(options.nclusters+1:end, 1:options.nclusters) = mu_val*eye(nstates/2);
        ss_A(options.nclusters+1:end, options.nclusters+1:end) = ...
            -1*mu_val*eye(nstates/2) + ss_A(options.nclusters+1:end, options.nclusters+1:end);
    end
    
    % state to observation matrix (assuming states are observations)
    ss_C = eye(nstates);
    % create process noise matrix 
    train_extrema = [train_maxlocs, train_minlocs];
    train_extrema = train_extrema(:, extrema_gmm_ord);
    train_mu = [reshape(pos_gmm.model.mu(:, 1), 1, []), reshape(neg_gmm.model.mu(:, 1), 1, [])];
    train_mu = train_mu(extrema_gmm_ord);
   
    if options.qmatrix_abs
        train_mu_vec = zeros(size(train_extrema));
        train_mu_vec(1, :) = train_mu;
        for i = 2:size(train_extrema, 1)
            train_mu_vec(i, :) = train_mu_vec(i-1, :) + (mu_val)*(train_extrema(i, :) - train_mu_vec(i-1, :));
        end
        %train_mu_vec = movmean(train_extrema, 10);
        covMatrix = abs(cov([train_extrema, train_mu_vec]));
    else

        covpeaks = diff(train_extrema);
        covmu = diff((train_extrema - train_mu)*mu_val);
        covMatrix = cov([covpeaks, covmu]);
        %covMatrix = diag(diag(covMatrix)); % edit 2/17/2022
    end
    ss_S = zeros(size(ss_A, 1), size(ss_C, 1)); 

    % set the initial state...but this will definitely be different 
    kFstatePred = [train_extrema(1, :), train_mu];
    
elseif options.singlekalman
    
    nstates = options.nclusters; % based on number of kalman filters prior ]  

    % predicted big boi
    P_pred = zeros(nstates, nstates);

    % residual big boi 
    res = zeros(nbeats, nstates);
    
    % state estimation big boi
    stateEst = zeros(nbeats, nstates);  
    
    % state transition matrix, state to obs matrix 
    ss_A = eye(nstates);
    ss_C = eye(nstates);
    ss_S = zeros(size(ss_A, 1), size(ss_C, 1));
    
    % 
    train_extrema = [train_maxlocs, train_minlocs];
    train_extrema = train_extrema(:, extrema_gmm_ord);
    covpeaks = diff(train_extrema);
    covMatrix = cov(covpeaks);
    %covMatrix = diag(diag(covMatrix));
    kFstatePred = train_extrema(1, :);
    
end
    
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   2a Prepping Data for Testing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% find specified features for each scg beat 
if options.ndim == 2
    [pos_cand, neg_cand, pos_prom, neg_prom] = findAO(scg_beats, options.ncand, fs,...
        [options.axrange(1), options.axrange(2)], 1);
elseif options.ndim == 1
    [pos_cand, neg_cand] = findAO(scg_beats, options.ncand, fs,...
        [options.axrange(1), options.axrange(2)], 1);
end

% set the features for each peak/valley
neg_data = zeros([nbeats, options.ncand, options.ndim]);
neg_data(:, :, 1) = neg_cand;
pos_data = zeros([nbeats, options.ncand, options.ndim]);
pos_data(:, :, 1) = pos_cand;

% if 2 extracted features add in prominence 
if options.ndim == 2
    pos_data(:, :, 2) = pos_prom;
    neg_data(:, :, 2) = neg_prom;
end

% for each testing beat 
for beat = 1:size(pos_cand, 1)
    if beat == 785
        beat;
    end
    
    % prep positive data for GMM
    [test_posdata, ~] = prep_data(pos_data, beat);
    [~, unique_idx] = unique(test_posdata(:, 1));
    test_posdata = test_posdata(unique_idx, :);
    
    % prep negative data for GMM
    [test_negdata, ~] = prep_data(neg_data, beat);
    [~, unique_idx] = unique(test_negdata(:, 1));
    test_negdata = test_negdata(unique_idx, :);
    
    % sort the peaks/valley and their features by their location
    [~, cand_sort_idx ] = sort([test_posdata(:, 1);test_negdata(:, 1)],'ascend');
    
    % if no data is available skip! holy.... if we skip (which we don't
    % really, but if we do it looks at previous beats which currently have
    % 0??)
    if sum(isnan(test_posdata)) + sum(isnan(test_negdata)) > 0
        fprintf(['Beat ', num2str(beat), ' has no features and will be skippe\n'])
        continue;
    end
    
    % evaluate peak and valley data on positive gmm 
    [~, pos_post_pdf, ~, pos_d2] = pos_gmm.evalPoints(test_posdata, 'type', 'mdist', 'delay', del);
    [~, ~, ~, pos_neg_d2] = pos_gmm.evalPoints(test_negdata, 'type', 'mdist', 'delay', del);
    
    % evaluate peak and valley data on negative gmm
    [~, neg_post_pdf, ~, neg_d2] = neg_gmm.evalPoints(test_negdata, 'type', 'mdist', 'delay', del);
    [~, ~, ~, neg_pos_d2] = neg_gmm.evalPoints(test_posdata, 'type', 'mdist', 'delay', del);
    
    % if using the relative distances between a peak/valley to all clusters
    %
    % break down the distnace in terms of time and prom dimensions and ONLY
    % alter time if using multiple dimensions 
    if options.relweights
        
        pos_d2_time = pos_d2;
        pos_neg_d2_time = pos_neg_d2;
        
        neg_d2_time = neg_d2;
        neg_pos_d2_time = neg_pos_d2;
        
        if options.ndim == 2
            % evaluate peak and valley data on positive gmm time dim
            [~, ~, ~, pos_d2_time] = pos_gmm.evalPoints(test_posdata, 'type', 'mdist', 'delay', 1);
            [~, ~, ~, pos_neg_d2_time] = pos_gmm.evalPoints(test_negdata, 'type', 'mdist', 'delay', 1);

            % evaluate peak and valley data on negative gmm time dim
            [~, ~, ~, neg_d2_time] = neg_gmm.evalPoints(test_negdata, 'type', 'mdist', 'delay', 1);
            [~, ~, ~, neg_pos_d2_time] = neg_gmm.evalPoints(test_posdata, 'type', 'mdist', 'delay', 1);
            
            % calculate distance in terms of prominence
            pos_d2_prom = sqrt(pos_d2.^2 - pos_d2_time.^2);
            pos_neg_d2_prom = sqrt(pos_neg_d2.^2 - pos_neg_d2_time.^2);
            
            %
            neg_d2_prom = sqrt(neg_d2.^2 - neg_d2_time.^2);
            neg_pos_d2_prom = sqrt(neg_pos_d2.^2 - neg_pos_d2_time.^2);
            
        end
        % placeholder for peaks and valley distances to all gmm cluster means
        pos_adj_d2 = zeros(size(pos_d2));
        neg_adj_d2 = zeros(size(neg_d2));
        
        % concatenate peak and valley distances to all gmm clusters and
        % sort by gmm cluster locations
        pos_cand_d2 = [pos_d2_time, neg_pos_d2_time]; pos_cand_d2 = pos_cand_d2(:, extrema_gmm_ord);
        neg_cand_d2 = [pos_neg_d2_time, neg_d2_time]; neg_cand_d2 = neg_cand_d2(:, extrema_gmm_ord);      
        
        % combine all distances together and sort by peak/valley locations
        all_gmm_d2 = [pos_cand_d2;neg_cand_d2];
        all_gmm_d2 = all_gmm_d2(cand_sort_idx, :);
        
        % for each peak cluster 
        for i = 1:npkclusters
            % take the mean of the abs differences of peak to all gmm clusters
            % to the ith gmm peak cluster to all other clusters
            pos_adj_d2(:, i) = mean(abs(pos_cand_d2 - unprep_pos_relations(:, i)'), 2);
        end
        
        % for each valley cluster
        for i = 1:nvalclusters
            % take the mean of the abs differences of valleys to all gmm clusters
            % to the ith gmm valley cluster to all other clusters
            neg_adj_d2(:, i) = mean(abs(neg_cand_d2 - unprep_neg_relations(:, i)'), 2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save original distances (for debugging)
        pos_d2_temp = pos_d2;
        neg_d2_temp = neg_d2;
        
        % overwrite with adjusted distances with weighted normalization
        pos_d2 = pos_adj_d2;
        neg_d2 = neg_adj_d2;
        
        % if using prom and time dimensions update the total distance with
        % the increased or reduced distances 
        if options.ndim == 2
            
            % new d2 is sqrt(
            pos_d2 = sqrt(pos_adj_d2.^2 + pos_d2_prom.^2);
            neg_d2 = sqrt(neg_adj_d2.^2 + neg_d2_prom.^2);
            
        end
       
    end
    
    % normalize p(g_i|x_j)/p(g_i|x)
    pos_gnpost_pdf = pos_post_pdf ./ sum(pos_post_pdf, 1);
    neg_gnpost_pdf = neg_post_pdf ./ sum(neg_post_pdf, 1);
     
    % get the total number of candidates and number of peak/valley 
    npkcands = size(test_posdata, 1);
    nvalcands = size(test_negdata, 1);
    nextrema = nvalcands+npkcands;

    % save data, d2, probs, nrobs into struct array info for easy viewing
    % later 
    test_exttime = [test_posdata(:, 1);test_negdata(:, 1)];
    test_extsign = [ones(npkcands, 1);zeros(nvalcands, 1)];         
    % sort the locations and apply learned sorting to signs 
    [~, extrema_cand_ord] = sort(test_exttime, 'ascend');
    test_extdata = [test_posdata; test_negdata];
    test_extdata = test_extdata(extrema_cand_ord, :);
    test_extsign = test_extsign(extrema_cand_ord);  
    % matrix to use for prob saving 
    prob_matrix = zeros(nextrema, options.nclusters);
    prob_matrix(1:npkcands, 1:npkclusters) = pos_post_pdf;
    prob_matrix(end-nvalcands+1:end, end-nvalclusters+1:end) = neg_post_pdf;
    prob_matrix = prob_matrix(extrema_cand_ord, :);
    prob_matrix = prob_matrix(:, extrema_gmm_ord);
    % matrix to use for prob saving 
    d2_matrix = ones(nextrema, options.nclusters)*inf;
    d2_matrix(1:npkcands, 1:npkclusters) = pos_d2;
    d2_matrix(end-nvalcands+1:end, end-nvalclusters+1:end) = neg_d2;
    d2_matrix = d2_matrix(extrema_cand_ord, :);
    d2_matrix = d2_matrix(:, extrema_gmm_ord);
    % matrix to use for prob saving 
    nprob_matrix = zeros(nextrema, options.nclusters);
    nprob_matrix(1:npkcands, 1:npkclusters) = pos_gnpost_pdf;
    nprob_matrix(end-nvalcands+1:end, end-nvalclusters+1:end) = neg_gnpost_pdf;
    nprob_matrix = nprob_matrix(extrema_cand_ord, :);
    nprob_matrix = nprob_matrix(:, extrema_gmm_ord);
    info(beat).d2 = d2_matrix;
    info(beat).probs = prob_matrix;
    info(beat).nprobs = nprob_matrix;
    info(beat).data = [test_extsign, test_extdata];
    
    
    
    % placeholder for sorted peaks, distances, and probabilities 
    pos_sortedPeaks_i = zeros(options.ndim, npkclusters);
    pos_sortedPeaks_d2_i = zeros(1, npkclusters);
    pos_sortedPeaks_prob_i = zeros(1, npkclusters);
    pos_sortedPeaks_nprob_i = zeros(1, npkclusters);
    
    % placeholder for sorted valleys, distances, and probabilities
    neg_sortedPeaks_i = zeros(options.ndim, nvalclusters);
    neg_sortedPeaks_d2_i = zeros(1, nvalclusters);
    neg_sortedPeaks_prob_i = zeros(1, nvalclusters);
    neg_sortedPeaks_nprob_i = zeros(1, nvalclusters);
    
    % for each peak cluster
    for i = 1:npkclusters
        
        % get the peak features that correpsond to clsuter i 
        cluster_i_cand = test_posdata;
        
        % get the associated distances/probs for peaks in cluster i 
        d2_i_cand = pos_d2(:, i)';
        prob_i_cand = pos_post_pdf(:, i)';
        nprob_i_cand = pos_gnpost_pdf(:, i)';
        
        % find the idx in the cluster values taht has the best parameters 
        cluster_i_idx_prob = pos_gmm.findCluster(d2_i_cand, 'min');
        %cluster_i_idx_prob = pos_gmm.findCluster(prob_i_cand, 'max');

        % save the peaks, distances, and probabilities
        pos_sortedPeaks_i(:, i) = cluster_i_cand(cluster_i_idx_prob, :);
        pos_sortedPeaks_d2_i(i) = d2_i_cand(cluster_i_idx_prob);
        pos_sortedPeaks_nprob_i(i) = nprob_i_cand(cluster_i_idx_prob);
        
    end 
    
    % for each valley cluster
    for i = 1:nvalclusters
        
        % get the valley features that correpsond to valley cluster i 
        cluster_i_cand = test_negdata;
        
        % get the associated distances/probs for valleys in cluster i 
        d2_i_cand = neg_d2(:, i)';
        prob_i_cand = neg_post_pdf(:, i)';
        nprob_i_cand = neg_gnpost_pdf(:, i)';
        
        % find the idx in the cluster values taht has the best parameters 
        cluster_i_idx_prob = pos_gmm.findCluster(d2_i_cand, 'min');
        %cluster_i_idx_prob = neg_gmm.findCluster(prob_i_cand, 'max');

        % save the valleys, distances, and probabilities 
        neg_sortedPeaks_i(:, i) = cluster_i_cand(cluster_i_idx_prob, :);
        neg_sortedPeaks_d2_i(i) = d2_i_cand(cluster_i_idx_prob);
        neg_sortedPeaks_prob_i(i) = prob_i_cand(cluster_i_idx_prob);
        neg_sortedPeaks_nprob_i(i) = nprob_i_cand(cluster_i_idx_prob);
        
    end
    
    % concatentate the peak and valley features, distances, and probabilities
    sortedPeaks_i = [pos_sortedPeaks_i(1, :), neg_sortedPeaks_i(1, :)];
    sortedPeaks_d2_i = [pos_sortedPeaks_d2_i, neg_sortedPeaks_d2_i];
    sortedPeaks_prob_i = [pos_sortedPeaks_prob_i, neg_sortedPeaks_prob_i];
    sortedPeaks_nprob_i = [pos_sortedPeaks_nprob_i, neg_sortedPeaks_nprob_i];
    
    % sort the concatenated features, distance,s and probs based on
    % location and save
    sortedPeaks(beat, :) = sortedPeaks_i(extrema_gmm_ord);
    sortedPeaks_d2(beat, :) = sortedPeaks_d2_i(extrema_gmm_ord);
    sortedPeaks_probs(beat, :) = sortedPeaks_prob_i(extrema_gmm_ord);
    sortedPeaks_nprobs(beat, :) = sortedPeaks_nprob_i(extrema_gmm_ord);
   
    if options.ndim == 2
        sortedProms_i = [pos_sortedPeaks_i(2, :), neg_sortedPeaks_i(2, :)];
        sortedProms(beat, :) = sortedProms_i(extrema_gmm_ord);
    end
    %%%%%%%%%%%%%%%%%%%%% SHOULD CHECK THIS FOR SOME REASON THE DISTANCES
    %%%%%%%%%%%%%%%%%%%%% ARE 0 SOMETIMES
    pos_d2(pos_d2 == 0) = 1e-3;
    neg_d2(neg_d2 == 0) = 1e-3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % run the Dijkkstra algorithm if specified 
    if options.dijkstra

        % get the total number of candidates and number of peak/valley 
        npkcands = size(test_posdata, 1);
        nvalcands = size(test_negdata, 1);
        nextrema = nvalcands+npkcands;
        
        % concatenate peak/valley locations and signs 
        test_extdata = [test_posdata(:, 1);test_negdata(:, 1)];
        test_extsign = [ones(npkcands, 1);zeros(nvalcands, 1)];         
        
        % sort the locations and apply learned sorting to signs 
        [test_extdata, extrema_cand_ord] = sort(test_extdata, 'ascend');
        test_extsign = test_extsign(extrema_cand_ord); 
        
        if sum(diff(test_extsign) == 0) > 0
            disp(['Beat: ', num2str(beat), ' Not alternating pk/valley'])
        end
        
        if options.ndim == 2
            prom_extdata = [test_posdata(:, 2); test_negdata(:, 2)];
            prom_extdata = prom_extdata(extrema_cand_ord);
        end
        
        % placeholder matrix based on normalized p(x|g) for all peak and valleys 
        test_matrix = zeros(nextrema, options.nclusters);
        
        % fill in pk to pk gmm and val to val gmm probabilities 
        test_matrix(1:npkcands, 1:npkclusters) = pos_gnpost_pdf;
        test_matrix(end-nvalcands+1:end, end-nvalclusters+1:end) = neg_gnpost_pdf;
        
        % sort based on candidate order and gmm order
        test_matrix = test_matrix(extrema_cand_ord, :);
        test_matrix = test_matrix(:, extrema_gmm_ord);
                
        % placeholder matrix based on distances for all pks and valleys
        % with inf distance for pk to val gmm and val to pk gmm
        value_matrix = ones(nextrema, options.nclusters)*inf; 
        
        if options.dijk_mixd2
            % fill in pk to pk gmm and val to val gmm probabilities 
            value_matrix(1:npkcands, 1:npkclusters) = arrayfun(@(x) min(x, 1e4), exp((1./pos_gnpost_pdf)-1).*pos_d2);
            value_matrix(end-nvalcands+1:end, end-nvalclusters+1:end) = arrayfun(@(x) min(x, 1e4), exp((1./neg_gnpost_pdf)-1).*neg_d2); 
        else
            % fill in pk to pk gmm and val to val gmm probabilities 
            value_matrix(1:npkcands, 1:npkclusters) = pos_d2;
            value_matrix(end-nvalcands+1:end, end-nvalclusters+1:end) = neg_d2; 
        end
        
        % sort based on caiddate and gmm order 
        value_matrix = value_matrix(extrema_cand_ord, :);
        value_matrix = value_matrix(:, extrema_gmm_ord);
        
        % remove certain nodes based on the distance from the means 
        if beat > 1
            mu_temp = [reshape(pos_gmm.mu(:, 1), 1, []), reshape(neg_gmm.mu(:, 1), 1, [])];
            mu_temp = mu_temp(extrema_gmm_ord);
%             mu_temp = sortedPeaks(beat-1, :);
%             extdata_temp = repmat(test_extdata, 1, 3);
%             extdata_temp = abs(extdata_temp - mu_temp);
%             extdata_temp_idx = (extdata_temp > 40*fs/1000);
%             value_matrix(extdata_temp_idx) = 1e3;
        end
        
        % matrix to use for prob saving 
        prob_matrix = zeros(nextrema, options.nclusters);
        prob_matrix(1:npkcands, 1:npkclusters) = pos_post_pdf;
        prob_matrix(end-nvalcands+1:end, end-nvalclusters+1:end) = neg_post_pdf;
        prob_matrix = prob_matrix(extrema_cand_ord, :);
        prob_matrix = prob_matrix(:, extrema_gmm_ord);
        
        % matrix to use for prob saving 
        nprob_matrix = zeros(nextrema, options.nclusters);
        nprob_matrix(1:npkcands, 1:npkclusters) = pos_gnpost_pdf;
        nprob_matrix(end-nvalcands+1:end, end-nvalclusters+1:end) = neg_gnpost_pdf;
        nprob_matrix = nprob_matrix(extrema_cand_ord, :);
        nprob_matrix = nprob_matrix(:, extrema_gmm_ord);
        
        
        % to acount for distances in the nodes considered in the first
        % clusters, add a column to the left and a row above to act as our "dummy node"
        test_matrix = [zeros(size(test_matrix, 1), 1), test_matrix];
        test_matrix = [zeros(1, size(test_matrix, 2)); test_matrix];
        test_matrix(1, 1) = 1;
        
        % add dummy node to distance matrix
        value_matrix = [zeros(size(value_matrix, 1), 1), value_matrix];
        value_matrix = [zeros(1, size(value_matrix, 2)); value_matrix];
        
        % create the dijkstra object and find the nodes 
        dijk = Dijkstra();
        dijk = dijk.setval('hweight', options.hweight);
        
        % add in dummy node that is the reverse extrema type of the first cluster
        dijk = dijk.setval('exttype', [~extrema_first;test_extsign]);
        dijk = dijk.prep_dijkstra(test_matrix, options.dijkthresh, value_matrix);

        % figure out all candidate node routes and corresponding path
        % distances'
        [node_dist, node_routes, node_dist_routes] = dijk.calc_paths(dijk.start_nodes, dijk.end_nodes);

        % need to alter column because the "dummy node" offsets the matrix by 1
        dijk = dijk.setval('colidx', dijk.colidx - 1);
        dijk = dijk.setval('rowidx', dijk.rowidx - 1);
        
        % convert node number to peak/valley index number and corresponding value
        pk_idx_routes = dijk.rowidx(node_routes(:, 2:end));
        pk_vals_routes = test_extdata(pk_idx_routes);
        
        % now find the best route using the min distance 
        [best_dist, best_route_idx] = min(node_dist);
        best_node_route = node_routes(best_route_idx, 2:end);
        
        % determine the best route if there are more than 1 route
        % candidates 
        custom_choice = 0;
        
        % if there are multiple candidate routes 
        if length(node_dist) ~= 1
            
            % testing region 
            if custom_choice && beat > 1
                % try 2 things:
                %   1. use previous peaks 
                %   2. if 2 
                % find top 3 routes (if available)
                [~, best_cand] = sort(node_dist, 'ascend');
                try
                % get the top 2 candidates routes based on distance 
                best_cand = best_cand(1:min(2, length(best_cand)));
                catch
                    best_cand;
                end
                
                % reduce the size of the routes in terms of 
                pk_vals_routes = pk_vals_routes(best_cand, :);
                pk_idx_routes = pk_idx_routes(best_cand, :);
                node_dist_routes = node_dist_routes(best_cand, :);
                
                % use prior sorted extremas as reference to choose a route
                [~, best_route_idx] = min(sum((sortedPeaks(beat - 1, :) - pk_vals_routes).^2, 2));
            end
            
            % get the location, indices, and distance of the pk/val sequence
            best_pk_vals_route = pk_vals_routes(best_route_idx, :);
            best_pk_idx_route = pk_idx_routes(best_route_idx, :);
            best_pk_dist_route = node_dist_routes(best_route_idx, :);
           
        % if there is only one route available 
        else
            best_pk_vals_route = pk_vals_routes';
            best_pk_idx_route = pk_idx_routes';
            best_pk_dist_route = node_dist_routes;
        end

        % separate into positive and negative peaks 
        best_pk_sign_route = test_extsign(best_pk_idx_route);
        
        % save the sorted peak estimations 
        sortedPeaks(beat, :) = best_pk_vals_route;
        if options.ndim == 2
            sortedProms(beat, :) = prom_extdata(best_pk_idx_route);
        end
        % if specified, save distances as dijkstra altered distances 
        if options.dijkweight
            sortedPeaks_d2(beat, :) = best_pk_dist_route(2:end);
            
        % otherwise use original distances 
        else
            sortedPeaks_d2(beat, :) = general.indexMatrix(value_matrix(2:end, 2:end), best_pk_idx_route')';
        end
        
        % save probabilities 
        sortedPeaks_probs(beat, :) = general.indexMatrix(prob_matrix, best_pk_idx_route')';
        sortedPeaks_nprobs(beat, :) = general.indexMatrix(nprob_matrix, best_pk_idx_route')';
        
    end
    
    
    if options.multikalman
        
        % here we are assuming the observation/state vector matrix is  
        % [p1, p2, p3, p4, u1, u2, u3, u4];
        % kF obs dependent on which peaks you are choosing 
        peakobs = sortedPeaks(beat, :);
        
        % extract and sort current location means of pk/val gmm 
        muobs = [reshape(pos_gmm.mu(:, 1), 1, []), reshape(neg_gmm.mu(:, 1), 1, [])];
        muobs = muobs(extrema_gmm_ord);
        
        % if we want to adjust the mean observation 
        if options.prealter
             muobs = muobs + (peakobs - muobs).*sortedPeaks_nprobs(beat, :);
        end
        
        % create state vector, combining extrema vector and means 
        kFobs = [peakobs, muobs];
        
        % given the best route in nodes, find the associated diatances 
        %node_ind = sub2ind(size(d2), dijk.rowidx(best_node_route), ...
        %    dijk.colidx(best_node_route))';
        
         % set the observation nosie covariacne R using the mahalbonis
         % distances and upper threshold if necessary 
         %peakR = d2(node_ind);
%         [~, ~, ~, pos_d2] = pos_gmm.evalPoints(test_posdata, 'type', 'mdist', 'delay', del);
%         [~, ~, ~, pos_neg_d2] = pos_gmm.evalPoints(test_negdata, 'type', 'mdist', 'delay', del);
%         [~, neg_post_pdf, ~, neg_d2] = neg_gmm.evalPoints(test_negdata, 'type', 'mdist', 'delay', del);
%         [~, ~, ~, neg_pos_d2] = neg_gmm.evalPoints(test_posdata, 'type', 'mdist', 'delay', del);

         % calculate distances JUST on time
         [~, ~, ~, pos_d2] = pos_gmm.evalPoints(peakobs(extrema_signs == 1), 'type', 'mdist', 'delay', 1);
         [~, ~, ~, neg_d2] = neg_gmm.evalPoints(peakobs(extrema_signs == 0), 'type', 'mdist', 'delay', 1);
         
         % use the distances saved 
         if options.mixd2
             peakR = exp(1./sortedPeaks_nprobs(beat, :) - 1).*sortedPeaks_d2(beat, :)*options.rconst;
         else
            peakR = sortedPeaks_d2(beat, :)*options.rconst;
         end
         % if using distances JUST on time replace R vector with those
         % distances
         %peakR(extrema_signs == 1) = pos_d2; peakR(extrema_signs == 0) = neg_d2;
         
         % calculate 
         muR = peakR*mu_val; 
         R_vec = [peakR, muR];
         kFthresh = 1e4*ones(1, nstates);
         ss_R = diag(min(abs([R_vec; kFthresh])));
         
         % set the process noise covariance 
         ss_Q = options.qconst*covMatrix;
         
         info(beat).R = ss_R;
         info(beat).P = P_pred;
         % actually stateEst is the output or the state to output matrix C
         % times kFstatePred the posterior state prediction
         % for the updating single step mu 
         [stateEst(beat, :), kFstatePred, P_pred, res(beat, :)] = ...
             kalmanFilterAlt(ss_A, ss_C, ss_Q, ss_R, ss_S, ...
             kFobs, 'forLoop', kFstatePred, P_pred);
         info(beat).res = res(beat, :);   
         info(1).Q = ss_Q;

         % if are trying to impose some 
         if options.mkfadaptive
             pos_gmm.mu = pos_gmm.mu + (sortedPeaks(beat, extrema_signs==1)' - pos_gmm.mu).*1./peakR(extrema_signs == 1)'*mu_val;
             neg_gmm.mu = neg_gmm.mu + (sortedPeaks(beat, extrema_signs==0)' - neg_gmm.mu).*1./peakR(extrema_signs == 0)'*mu_val;
             
         % if     
         else
             
             % extract new means from state equations and replace pk/val
             % gmm location means 
             new_mu = stateEst(beat, options.nclusters+1:end)';
             pos_gmm.mu(:, 1) = new_mu(extrema_signs == 1);
             neg_gmm.mu(:, 1) = new_mu(extrema_signs == 0);
             
         end
         if options.pkadj
             sortedPeaks(beat, :) = stateEst(beat, 1:options.nclusters);
         end

    end
    
    if options.singlekalman 
    
        % kF obs dependent on which peaks you are choosing 
        kFobs = sortedPeaks(beat, :);

        % given the best route in nodes, find the associated diatances 
        % node_ind = sub2ind(size(d2), dijk.rowidx(best_node_route), ...
        %    dijk.colidx(best_node_route))';

         % set the observation nosie covariacne R using the mahalbonis
         % distances and upper threshold if necessary 
         % R_vec = d2(node_ind);
         % use the distances saved 
         if options.mixd2
             peakR = exp(1./sortedPeaks_nprobs(beat, :) - 1).*sortedPeaks_d2(beat, :)*options.rconst;
         else
            peakR = sortedPeaks_d2(beat, :)*options.rconst;
         end
         R_vec = peakR;
         kFthresh = 1e4*ones(1, nstates);
         ss_R = diag(min(abs([R_vec; kFthresh])));
         
         % set the process noise covariance matrix 
         ss_Q = options.qconst*covMatrix;
         
         % for the 
         [stateEst(beat, :), kFstatePred, P_pred, res(beat, :)] = ...
             kalmanFilterAlt(ss_A, ss_C, ss_Q, ss_R, ss_S, ...
             kFobs, 'forLoop', kFstatePred, P_pred);

          sortedPeaks(beat, :) = stateEst(beat, 1:options.nclusters);
         % update the distances with the new peaks provided 
         %%%%% 10/25/2021 Since this is smoothening i don't actually want
         %%%%% to do this because I want the distance of the original peaks
         %%%%% 
         %[~, newpost_pdf, ~, newd2] = gmm.evalPoints(sortedPeaks(beat, :)', 'type', 'mdist');
         %sortedPeaks_d2(beat, :) = newd2';
         %sortedPeaks_prob(beat, :) = newpost_pdf'; 
        
    end
    
    % update 
    if options.gmmupdate && ~options.multikalman
        
        % update mu 
        old_pos_mu = pos_gmm.mu(:, 1);
        
        new_pks = sortedPeaks(beat, extrema_signs == 1)';
        pos_gmm.mu(:, 1) = old_pos_mu + (new_pks - old_pos_mu) * mu_val; 
        
        old_neg_mu = neg_gmm.mu(:, 1);
        new_vals = sortedPeaks(beat, extrema_signs == 0)';
        neg_gmm.mu(:, 1) = old_neg_mu + (new_vals - old_neg_mu) * mu_val; 
            
    end
    muobs = [reshape(pos_gmm.mu(:, 1), 1, []), reshape(neg_gmm.mu(:, 1), 1, [])];
    muobs = muobs(extrema_gmm_ord);
    allmu(beat, :) = muobs;
         
end

if options.plot
    figure;
    ax(1) = subplot(2, 1, 1);
    imagesc(scg_beats); newmap = contrast(scg_beats); colormap(newmap); hold on
    for i = 1:(options.nclusters)
        p = plot(sortedPeaks(:, i), 'LineWidth', 2);
        %p.Color(4) = 0.5;
    end
    
    xlabel('Beat Number')
    ylabel('Delay (samples)')
    ax(2) = subplot(2 ,1 ,2);
    plot(sortedPeaks_nprobs); title('Probabilities (against other candidates)')
    linkaxes(ax, 'x')

end

% at the end merge proms with peaks 
if options.ndim == 2
    sortedPeaks = cat(3, sortedPeaks, sortedProms);
end
varargout{1} = allmu;
varargout{2} = sortedPeaks_d2;
varargout{3} = sortedPeaks_probs;
varargout{4} = sortedPeaks_nprobs;
varargout{5} = info;
%varargout{2} = pos_post_pdf;
%varargout{3}= pos_gmm;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Helper Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_locs, best_locs_sign] = find_extremas(template, ax_range, fs, perc_lim, Plot)
% find the best pk/valley sequence to end up tracking using 
    % get the peaks and valley locations, amplitudes, and prominences 
    
    ax_range = ax_range*fs/1000;
    ax_range = ax_range(1):ax_range(end);
    [ppks, plocs, ~, pproms] = findpeaks(template(ax_range), 'Annotate','extent');
    [~, nlocs, ~, nproms] = findpeaks(template(ax_range)*-1);
    
    % adjust locations to beginning of signal
    plocs = plocs + ax_range(1);
    nlocs = nlocs + ax_range(1);
    npks = template(nlocs);


    % for plotting, plot the pks/valleys with highest prominences 
    [~, pproms_ord] = sort(pproms, 'descend');
    [~, nproms_ord] = sort(nproms, 'descend');
    
    red = [1, 0, 0];
    blue = [0, 0, 1];
    nlim = 3;
    if Plot
        figure; plot(template); hold on;
        scatter(plocs(pproms_ord(1:3)), template(plocs(pproms_ord(1:3))), [], [repmat(red, nlim, 1)'.*linspace(1, 0.4, nlim)]', 'filled');
        scatter(nlocs(nproms_ord(1:3)), template(nlocs(nproms_ord(1:3))), [], [repmat(blue, nlim, 1)'.*linspace(1, 0.4, nlim)]', 'filled');
        title('template with top 3 peaks/valleys labeled based on prominence')
        legend('template', 'peaks', 'valleys')
    end
    % combine pk/valley locations and sort based on location
    locs_seq = [plocs; nlocs];
    [locs_seq, locs_ord] = sort(locs_seq, 'ascend');
    
    % combined prominences, pk sign, amplitude based on prior location sort
    prom_seq = [pproms; nproms]; prom_seq = prom_seq(locs_ord);
    sign_seq = [ones(size(plocs)); zeros(size(nlocs))];  sign_seq = sign_seq(locs_ord);
    pks_seq = [ppks; npks]; pks_seq = pks_seq(locs_ord);
    
    % get mean of the amplitude difference to direct adjacent extrema
    avg_amp_seq = zeros(size(pks_seq));
    avg_amp_seq(2:end-1) = movmean(abs(diff(pks_seq)), [0, 1], 'Endpoints', 'discard')*2;
    
    % check if the sorted sequence alternates between peaks and valleys
    if sum(diff(sign_seq) == 0) ~= 0
        error('not alternating pos and neg peaks');
    end
    
    % checking sequences from 2 to 5 peak/valleys
    seq_lens = 2:5;
    
    % prominence placeholder and sequence 
    best_prom = zeros(length(seq_lens), 1);
    best_seq = cell(length(seq_lens), 1);
    
    % for each sequence length
    for seq_len = seq_lens
        
        % calculate the average prominence of adjacent peak/valley sequence
        prom_avg = movmean(avg_amp_seq, [0, seq_len-1], 'Endpoints', 'discard');
        
        % find the window that has the highest promince
        best_prom_idx = find(prom_avg == max(prom_avg));
        
        % save the prominence and sequence 
        best_prom(seq_len-1) = prom_avg(best_prom_idx);
        best_seq{seq_len-1} = best_prom_idx:best_prom_idx+seq_len-1;
    end

    % calculate percent reduction from highest prominence 
    perc_change = ((best_prom - best_prom(1))/best_prom(1))'*100;

    
    % find longest sequence with a reduction less than limit of largest
    % prominence
    best_seq_idx = find(abs(perc_change) < perc_lim*100, 1, 'last');
    
    % save the sequence indices 
    best_seq = best_seq{best_seq_idx};
    
    % output the locations 
    best_locs = locs_seq(best_seq);
    best_locs_sign = sign_seq(best_seq);
end

function [datavec, groupvec, org_dim] = prep_data(data, range)
    % here we assume data is a N x ncluster x ndim dataset 
    % makes a group vector Mxndim where M <= N 
    % makes a data vector Mxndim  where M <= N
    % range specifies the region to concatenate ex: 1:300
    
    
    % get the range of data specified
    data = data(range, :, :);
    
    % get the dimensions 
    [N, nclusters, ndims] = size(data);
    
    % allocate data vector with collapsed cluster dimension
    datavec = zeros(N*nclusters, ndims);

    % create group identifier with collapsed cluster dimension (for
    % training)
    group = repmat(1:nclusters, N, 1);
    groupvec = reshape(group, [], 1);

    % for each dimension collapse the data into Nx1 array
    for i = 1:ndims
        datavec(:, i) = reshape(data(:, :, i), [], 1);
    end

    % save the original dimensions of the data
    org_dim = [N, nclusters, ndims];
end

function [pos_cand, varargout] = findAO(scgbeats, num_cand, fs, ao_range, prom_flg)
    %%%%%%%%%%%%%% Calculate top 10 peaks in each beat within 250ms %%%%%%%%%%
    
    % initialize locations, prominences, and widths of relevants peaks and
    % valleys 
    pos_cand = zeros(size(scgbeats, 2), num_cand);
    neg_cand = zeros(size(scgbeats, 2), num_cand);
    neg_proms = neg_cand;
    pos_proms = pos_cand;
    pos_wids = pos_cand;
    neg_wids = neg_cand;
    
    % dictate a range to find peaks and valleys
    ao_range = ao_range*fs/1000;
    
    % for each beat
    for beat = 1:size(scgbeats, 2)
        
        % isolate the beat
        temp_beat = scgbeats(:, beat);
        
        % find prominence if desired 
        if prom_flg
            [pos_pks, pos_locs, pos_wids_i, pos_proms_i] = ...
                findpeaks(temp_beat(ao_range(1):ao_range(2)), 'MinPeakDistance', 0);%10*fs/1000);
            
        % otherwise find just locations 
        else
            [pos_pks, pos_locs] = findpeaks(temp_beat(ao_range(1):ao_range(2)), 'MinPeakDistance', 0);%10*fs/1000);
        end
        
        % if we manage to find peaks 
        if ~isempty(pos_pks)
            
            % realign location to begining of beat 
            pos_locs = pos_locs + ao_range(1) - 1;
            
            % if number of peaks less than number desired 
            if length(pos_locs) < num_cand
                
                % pad locations and peaks to number desired using the last
                % value
                pos_locs = [pos_locs; ones(num_cand - length(pos_locs), 1)*pos_locs(end)];
                pos_pks = [pos_pks; ones(num_cand - length(pos_pks), 1)*pos_pks(end)];
                
                % pad prominences and peaks if prominence specified
                if prom_flg
                    pos_proms_i = [pos_proms_i; ones(num_cand - length(pos_proms_i), 1)*pos_proms_i(end)];
                    pos_wids_i = [pos_wids_i; ones(num_cand - length(pos_wids_i), 1)*pos_wids_i(end)];
                end
                
            end
            
            % sort the peaks by prominence or actual peak valley
            %[~, pos_idx] = sort(pos_prom, 'descend');
            [~, pos_idx] = sort(pos_pks, 'descend');
            
            % get the top cand peaks and sort the top locations by location
            [pos_cand(beat, :), idx] = sort(pos_locs(pos_idx(1:num_cand)), 'ascend');
            
            % sort the prominences and widths if specified
            if prom_flg
                
                pos_proms_i = pos_proms_i(pos_idx(1:num_cand));
                pos_proms(beat, :) = pos_proms_i(idx);
                
                pos_wids_i = pos_wids_i(pos_idx(1:num_cand));
                pos_wids(beat, :) = pos_wids_i(idx);
                                
            end
        
        % if no peaks were found 
        else
            
            % fill the candidates with nans 
            pos_cand(beat, :) = nan*zeros(num_cand, 1)';
            
        end
         
        % find prominence and widths of valleys if speicifed 
        if prom_flg
            [~, neg_locs, neg_wids_i, neg_proms_i] = ...
                findpeaks(-1*temp_beat(ao_range(1):ao_range(2)), 'MinPeakDistance', 0);%10*fs/1000);
        % otherwise just find locations 
        else
            [~, neg_locs] = findpeaks(-1*temp_beat(ao_range(1):ao_range(2)), 'MinPeakDistance', 0);%10*fs/1000);
        end
        
        % if valleys were detected 
        if ~isempty(neg_locs)
            
            % realign locations to beginning of beat 
            neg_locs = neg_locs + ao_range(1) - 1;
            neg_pks = temp_beat(neg_locs);
            
            % if number of valleys less than desired number 
            if length(neg_locs) < num_cand
                
                % pad locations and peaks with last value
                neg_locs = [neg_locs; ones(num_cand - length(neg_locs), 1)*neg_locs(end)];
                neg_pks = [neg_pks; ones(num_cand - length(neg_pks), 1)*neg_pks(end)];
                
                % pad widths and prominences if specified 
                if prom_flg
                    neg_proms_i = [neg_proms_i; ones(num_cand - length(neg_proms_i), 1)*neg_proms_i(end)];
                    neg_wids_i = [neg_wids_i; ones(num_cand - length(neg_wids_i), 1)*neg_wids_i(end)];
                end
            end
            
            % sort valleys by prominence or value 
            %[~, neg_idx] = sort(neg_prom, 'descend');
            [~, neg_idx] = sort(neg_pks, 'ascend');
            
            % sort the locations of the top valleys 
            [neg_cand(beat, :), idx] = sort(neg_locs(neg_idx(1:num_cand)), 'ascend');
            
            % use the location sorting idx to sort prominence and widths 
            if prom_flg
                neg_proms_i = neg_proms_i(neg_idx(1:num_cand));
                neg_proms(beat, :) = neg_proms_i(idx);
                
                neg_wids_i = neg_wids_i(neg_idx(1:num_cand));
                neg_wids(beat, :) = neg_wids_i(idx);
            end
        
        % if no peaks were found 
        else
            % fill locations with nans 
            neg_cand(beat, :) = nan*zeros(num_cand, 1)';
        end

        

    end
    
    % output candidates, prominences, and widths 
    varargout{1} = neg_cand;
    varargout{2} = pos_proms;
    varargout{3} = neg_proms;
    varargout{4} = pos_wids;
    varargout{5} = neg_wids;
end

function [x_hat, varargout] = kalmanFilterAlt(A, C, Q, R, S, y, varargin)
% This function takes system matrices A, C and covariance matrices Q, R, S,
% and performs Kalman filtering on the observations y, to produce state
% estimates at each corresponding time point, x_hat
% This function can be used for as many steps as desired, as it simply
% iterates through the rows present in y

% Inputs
% A: state transition matrix - N x N
% C: state to output mapping matrix - M x N
% Q: process noise covariance - N x N
% R: measurement noise covariance - M x M
% S: cross-covariance between two noises - N x M
% y: output sequence - T x M, where T is number of time instances
% (optional)
% 'forLoop' - flag to specify whether for loop approach is being used
% x_pred_in: predicted state from previous step
% P_pred_in: predicted error covariance from previous step

% Outputs
% x_hat: estimated state sequence - T x N
% (optional)
% x_pred_out: predicted state for next step (varargout{1})
% P_pred_out: predicted error covariance for next step (varargout{2})

% Initialize default
forLoop = false;
update = 0; 

% Parse varargin
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'forLoop')
            forLoop = true;
            x_pred_in = varargin{arg + 1};
            P_pred_in = varargin{arg + 2};
        end
        if strcmp(varargin{arg}, 'Update')
            update = 1;
        end
    end
end

% System Equations (for reference)
% x_i+1 = A*x_i + w_i
% y_i = C*x_i + v_i
% E[w_i*w_i'] = Q
% E[v_i*v_i'] = R
% E[w_i*v_i'] = S

% Herein, we assume S = 0
if S(1) ~= 0
    disp('Function was not made for correlated noises')
end


% Initialize prediction and measurement update state vectors
x_hat = zeros(size(y, 1), size(A, 1));
x_pred = zeros(size(y, 1), size(A, 1));

% Initialize error covariance matrix
P_hat = zeros(size(y, 1), size(A, 1), size(A, 1));
P_pred = zeros(size(y, 1), size(A, 1), size(A, 1));


if forLoop
    % Initialize first prediction to the inputs to the function
    x_pred(1, :) = x_pred_in;
    P_pred(1, :, :) = P_pred_in;
else
    % Initialize first prediction to first observation, just to initialize
    % (essentially just implement naive prediction for first time step)
    x_pred(1, :) = y(1, :);
    P_pred(1, :, :) = zeros(size(A, 1), size(A, 1));
end


% For every timestep provided
for t = 1:size(y, 1)
    
    
    % Prediction update
    x_pred(t, :) = (A*x_pred(t, :)')';
    P_pred(t, :, :) = A*squeeze(P_pred(t, :, :))*A' + Q;
    
    % Measurement update
    % (note that transposes had to be used in certain places to change the
    % row instances into column vectors)
    Kalm_gain = squeeze(P_pred(t, :, :))*C'*inv(C*squeeze(P_pred(t, :, :))*C' + R);
    x_pred(t, :) = (x_pred(t, :)' + Kalm_gain*(y(t, :)' - C*x_pred(t, :)'))';
    P_pred(t, :, :) = (eye(size(A, 1)) - Kalm_gain*C)*squeeze(P_pred(t, :, :));
    x_hat(t, :) = C*x_pred(t, :)';

end

% The final prediction step actually appends to x_pred and P_pred beyond
% their original initializations. Hence, we can use that to our advantage
% if we would like to store that final prediction
if forLoop
    x_pred_out = x_pred(end, :);
    P_pred_out = squeeze(P_pred(end, :, :));
    residual = y(t, :)' - C*x_pred(t, :)';
    varargout{1} = x_pred_out;
    varargout{2} = P_pred_out; 
    varargout{3} = residual;
end
end

% Script but only works on David's computer rn 

folder = 'C:\Users\david\Dropbox (GaTech)\Saved Variables\Kalman';
pc = 'C:\Users\david\Dropbox (GaTech)';
irl = '/home/dlin/Dropbox (GaTech)';
%https://www.sciencedirect.com/science/article/pii/S0167865505003594#fig1
%%

type = irl;
filePath = [type, filesep, 'Saved Variables'];
savePath = [type, filesep, 'Dijtika_Plots_full_cov'];
seed = 420;
rng(seed);
%%
% Code outline 
% 1. Grab true RAO values  
% 2. Create noisy segments 
% 3. Train

% Two methods 
% just 1 observation (closest to previous estimate?) 

% take top 3 peaks send through 3 separate kalman filters 

sub = 1; 

for sub = 6%[1, 2, 3, 4, 5, 6]

%%%%% Create varying SNR SCG beat vector 
des_beats = 13000;
% pattern length 750 beats 
synpat = [10*ones(100, 1); 5*ones(50, 1); 0*ones(50, 1); -5*ones(50, 1); ...
    -5*ones(50, 1); -10*ones(50, 1); -15*ones(50, 1); -20*ones(100, 1); ...
    -15*ones(50, 1); -10*ones(50, 1); -5*ones(50, 1);...
    0*ones(50, 1); 5*ones(50, 1)];



rr_len = load([filePath, filesep, 'pig_data', num2str(sub), '.mat'], ...
    'rr');
rr_len = rr_len.rr;
des_beats = length(rr_len)-1;

rdict = repmat(synpat, ceil(des_beats/length(synpat)), 1);
rdict = rdict(1:des_beats);
if sub == 4
sepbeats_syn = zeros(1200, size(rdict, 1));
else
sepbeats_syn = zeros(1000, size(rdict, 1));
end
for snr = -20:5:10
%snr = 5;
disp(['Processing Subject: ', num2str(sub), ' at SNR: ', num2str(snr), 'dB'])
%%%%%%%%%%%%%%%%%%%%%% 1. Grab true RAO values  %%%%%%%%%%%%%%%%%%%%%%%



% first design a function to add noise 
scg = load([filePath, filesep, 'pig_data', num2str(sub), '.mat'], ...
    'scgz', 'Hd', 'rr', 'rao', 'rac');

% get the rr intervals, filter, rao, rac, scg (z axis)
rr = scg.rr;
Hd = scg.Hd;
rao = scg.rao;
rac = scg.rac;
scg = scg.scgz;

fs = 2000;
%%%%%%%%%%%%%%%%%%%%%% 2. Create noisy segments %%%%%%%%%%%%%%%%%%%%%%%

% parameters 
if sub == 4
    beat_len = 1200;
else
    beat_len = 1000;
end
noise_type = 'gauss';

% produce noise with specified snr
noise = prod_noise(scg, noise_type, seed);
adj_noise = noises.general.snr_noise(snr, scg, noise, 'variance');

% create noise signal 
noisy_scg = scg + adj_noise;

% preprocessing/separate into beats 
noisy_scg = filtfilt(Hd.Numerator, 1, noisy_scg);
sepbeats = cardio.general.separateBeats(noisy_scg, 'indices', rr, 'samples', beat_len, 'nanpad', 1);
sepbeats = sepbeats(:, 1:end-1);
%sepbeats = normalize(sepbeats, 1);

if snr == 10
    scg = filtfilt(Hd.Numerator ,1 , scg);
    sepbeats_cln = cardio.general.separateBeats(scg, 'indices', rr, 'samples', beat_len);
    sepbeats_cln = sepbeats_cln(:, 1:end-1);
    %sepbeats_cln = normalize(sepbeats_cln, 1);
    template = sum(sepbeats(:, 1:100), 2)/100;
end
sepbeats_syn(:, rdict == snr) = sepbeats(:, rdict == snr);
% create a template 


%figure(1);
%plot(template); title(['Template for Subject ', num2str(sub)])
%%


% 3 options population (will do this later), single template, sliding template
% do it via population vs sliding (for accurate sqi just take all sqi
% values above 0.4?
% just to load an sqi file 
sqi_file = load([filePath, filesep, 'Manifold', filesep, 'sqi_method_dtfm_noise_' , ...
    noise_type, '_snr' num2str(snr), '.mat'], 'sqi_tot');
sqi = sqi_file.sqi_tot{sub};

sqi_syn(rdict == snr) = sqi(rdict == snr);
end

% save original and normalized versions of the beats 
sepbeats = sepbeats_syn; 
unsepbeats = sepbeats; sepbeats = normalize(sepbeats, 1);
unsepbeats_cln = sepbeats_cln; sepbeats_cln = normalize(sepbeats_cln, 1);

% because for subject 1 the subject has 
if sub == 1
    unsepbeats = unsepbeats*-1;
    unsepbeats_cln = unsepbeats_cln * -1;
    sepbeats =sepbeats*-1;
    sepbeats_cln = sepbeats_cln *-1;
end
disp('Using sliding template to get features')

% here we find the "clean" (hopefully) tracking of the top 4 peaks using
% the clean signal 
[maxlocs, minlocs] = general.slideTemplate(sepbeats_cln, ...
    30, 30, fs, 'Closest', 'Axrange', [1, 250], ...
    'NumFeatures', 7, 'Emd', 'Template', mean(sepbeats_cln(:, 1:60), 2));
% here we find the tracking of the top 4 peaks using the noisy signal
[maxlocs_noisy, minlocs_noisy] = general.slideTemplate(sepbeats, ...
    30, 30, fs, 'Closest', 'Axrange', [1, 250], ...
    'NumFeatures', 7, 'Emd', 'Template', mean(sepbeats_cln(:, 1:60), 2));

% next we get we sort the locs by delay from r peaks for each signal 
maxlocs = maxlocs';
[~, sortidx] = sort(maxlocs(1, :), 'ascend');
maxlocs = maxlocs(:, sortidx);

% save the 
maxlocs_noisy = maxlocs_noisy';
[~, sortidx] = sort(maxlocs_noisy(1, :), 'ascend');
maxlocs_noisy = maxlocs_noisy(:, sortidx);


%%


% here we just plot what the peaks look like 
figure;hold on;
brange = 1:2000;
rao_amps = noises.general.indexMatrix(sepbeats_cln(:, brange), maxlocs(brange, :));
plot(sepbeats_cln(:, brange)); hold on;
for i = 1:size(rao_amps, 2)
    scatter(maxlocs(brange, i), rao_amps(brange, i))
end
title('Features found for clean signal')
figure;
plot(maxlocs)
%%

% trying number of dimensions: 
% ndim = 1 -> just rpk - peak time delay 
% ndim = 2 -> add in amplitude at that peak (TO DO: potentially try from
%                   the lower peak too 
% ndim = 3 -> add in peak width (TO DO: haven't implemented this yet) 
ndim = 1;

% number of peaks to track (with 7 tracked features this means 4 peaks)
nclusters = size(maxlocs, 2);

% create data frame with third dim just being each component of the signal 
data = zeros([size(maxlocs), ndim]);
data(:, :, 1) = maxlocs;
if ndim == 2
   data(:, :, 2) = rao_amps;
end

% initialize a gmm object using the cluster and dimensions 
gmm = GMM('nclusters', nclusters, 'ndim', ndim);

% estparams = 1;
% 
% % take out a subset of the data to act as the training set 
% train_range = 1:100;
% [train_data, train_labels] = prep_data(data, train_range);
% 
% % using the labels to enable a warm start for setting the iniital
% % means/sigmas for training 
% gmm = gmm.setStart(train_data, train_labels, estparams);
% gmm = gmm.fitModel(train_data);
% 
% %%
% 
% % take out a subset of the data to act as our testing set 
% test_range = 1:800;
% [test_data, test_labels, test_dims] = prep_data(data, test_range);
% 
% % evaluate the testing set 
% [idx, post_pdf, npost_pdf, d2] = gmm.evalPoints(test_data, 'type', 'gprob');
% 
% 
% % plot the differences in the predicted and true peak labels 
% figure; plot(idx); hold on; plot(test_labels)
% title('Small test set evaluation using GMM')
% legend('Predicted class', 'True Classs')
% 
% % plot the found distributions of the gmm with the data plotted as well
% wrong_labels = find((test_labels - idx) ~= 0); 
% gmm.genGMMcontour(test_data, 'points', test_data(wrong_labels, :)); 
% title('Trained GMM on Testing Data')
% xlabel('Delay'); ylabel('P(g_i, x_i)')

%
% %%%%%%%
% % 
% % Testing to see how the mean and covariances of the gmm changes through time 
% % 
% %%%%%%
% 
% % holds the models 
% models = {};
% 
% % holds the training/testing signals 
% signals = {};
% 
% % signal array length
% arr_len = 100;
% 
% % how many signals to offset the array 
% offset = 50;
% 
% % for each signal arry 
% for i = 1:offset:1800
%     
%     % get a subset of the data 
%     [train_data, train_labels] = prep_data(data, i:i+arr_len);
%     
%     % train a gmm model using the given points 
%     gmm = setStart(gmm, train_data, train_labels, estparams);
%     gmm = fitModel(gmm, train_data);
%     
%     % save the models and signals 
%     models = [models, {gmm}];
%     signals = [signals, {train_data}];
%     
% end
% % show a movie of the pdfs changing through time 
% gmmMovie(signals, models, 'speed', 'fast')
%

%%

% trying number of dimensions: 
% ndim = 1 -> just rpk - peak time delay 
% ndim = 2 -> add in amplitude at that peak (TO DO: potentially try from
%                   the lower peak too 
% ndim = 3 -> add in peak width / prominence(TO DO: haven't implemented this yet) 
ndim = 1;

% number of peaks to track (with 7 tracked features this means 4 peaks)
nclusters = size(maxlocs, 2);

% create data frame with third dim just being each component of the signal 
data = zeros([size(maxlocs), ndim]);
data(:, :, 1) = maxlocs;
if ndim == 2
   data(:, :, 2) = rao_amps;
end


% % initialize a gmm object using the cluster and dimensions 
gmm = GMM('nclusters', nclusters, 'ndim', ndim);


% take out a subset of the data to act as the training set 
train_range = 1:100;
[train_data, train_labels] = prep_data(data, train_range);

% using the labels to enable a warm start for setting the iniital
% means/sigmas for training 
gmm = gmm.setStart(train_data, train_labels);
gmm = gmm.fitModel(train_data);


%%%%%%%
% 
% Testing if for a signal we feed in the top 10 peaks can we still get the right AO tracked peaks  
% 
%%%%%%
% number of peaks to find for each signal 
num_cand = 10;
testclean = 1;
clean_vec = [1, 0];
for testclean = clean_vec
% determine if we are testing on 
if testclean
    disp('testing on clean beats')
    testing_beats = sepbeats_cln;
    untesting_beats = unsepbeats_cln;
    altmethod = maxlocs;
    cln_str = 'cln';
else
    disp('testing on noisy beats')
    testing_beats = sepbeats;
    untesting_beats = unsepbeats;
    altmethod = maxlocs_noisy;
    cln_str = 'noisy';
end


% candidates for each of the signals (first use the clean one to test)
rao_cand = findAO(testing_beats, num_cand, fs);

% find the associated amplitude for each candidate 
rao_amps = noises.general.indexMatrix(untesting_beats, rao_cand);

% create data structure for the GMM to use 
ao_data = zeros([size(rao_cand), ndim]);
ao_data(:, :, 1) = rao_cand;
if ndim > 1
    ao_data(:, :, 2) = rao_amps;
end

% hold solutions for different tests
dijknoupdate = zeros(size(rao_cand, 1), nclusters);
dijkupdate = zeros(size(rao_cand, 1), nclusters);
pkkfupdatekf = zeros(size(rao_cand, 1), nclusters); % kf smoothened pks for updates
pkkfupdatedi = zeros(size(rao_cand, 1), nclusters); % dijk pks for updates

% it goes [update; single kf; multikf]
tests = [0, 1, 1. 1;
            0, 0, 1, 0;
            0, 0, 0, 1];

for test_i = tests
% update parameter
update = test_i(1);
kFpeaks = test_i(2);
kFmulti = test_i(3);



% 4 scenarios 
if ~update && ~kFpeaks && ~kFmulti
    disp('using dijkojic with no updates and no kf')
elseif update && ~kFpeaks && ~kFmulti
    disp('using dijkojic with single step heuristic updates and no peak kf')
elseif update && kFpeaks && ~kFmulti
    disp('using dijkojic with single step heuristic updates using peak kf')
elseif update && ~kFpeaks && kFmulti
    disp('using dijkojic and peak and mu kalman updates')
else
    error('invalid combination')
end






% big boi stats 
% initialize the peaks 
%
nstates = nclusters; % based on number of kalman filters prior ]  

q = 0.2;
if kFpeaks
    
    nstates = nclusters; % based on number of kalman filters prior ]  

    % predicted big boi
    P_pred = zeros(nstates, nstates);

    % residual big boi 
    res = zeros(size(rao_cand, 1), nstates);

    saveobs = zeros(size(rao_cand, 1), nstates);
    % state estimation big boi
    stateEst = zeros(size(rao_cand, 1), nstates);    
    ss_A = eye(nstates);
    ss_C = eye(nstates);
    ss_S = zeros(size(ss_A, 1), size(ss_C, 1));
    covpeaks = diff(maxlocs(1:100, :));
    covMatrix = cov(covpeaks);
    %covMatrix = diag(diag(covMatrix));
    kFstatePred = maxlocs(1, :);
    
end

if kFmulti
    % here we assume state_vector = [p1, p2, p3, p4, u1, u2, u3, u4]
    nstates = nclusters*2; % based on number of kalman filters prior ]  
    mu_val = 1/100;
    % predicted big boi
    P_pred = zeros(nstates, nstates);

    % residual big boi 
    res = zeros(size(rao_cand, 1), nstates);

    saveobs = zeros(size(rao_cand, 1), nstates);
    % state estimation big boi
    stateEst = zeros(size(rao_cand, 1), nstates);    
    ss_A = eye(nstates);
    ss_A(5:end, 1:4) = mu_val*eye(nstates/2);
    ss_A(5:end, 5:end) = -1*mu_val*eye(nstates/2) + ss_A(5:end, 5:end);
    ss_C = eye(nstates);
     covpeaks = diff(maxlocs(1:100, :));
     covmu = diff((maxlocs(1:100, :) - reshape(gmm.model.mu, 1, []))*mu_val);
    covMatrix = cov([covpeaks, covmu]);
    ss_S = zeros(size(ss_A, 1), size(ss_C, 1));    
    kFstatePred = [maxlocs(1, :), reshape(gmm.model.mu, 1, [])];
    
    
end

% save 
all_dist = zeros(size(rao_cand, 1), nstates);
all_dist_multi = zeros(size(rao_cand, 1), nstates);

% create data structure that will hold the chosen peaaks among the
% candidates 
N = size(rao_cand, 1);
est_rao_dist = zeros(size(rao_cand, 1), nclusters);

est_rao_prob = est_rao_dist;
est_rao_mixed = est_rao_dist;
est_rao_dijk = est_rao_dist;
est_rao_multi = est_rao_dist;

% dummy to plot the output of one of the steps 
bstop = 1127;


% for each signal save the 
all_mu = est_rao_dist;
% brange = 2 already messing up?? 



for brange = 1:size(rao_cand, 1)
    
    if brange == bstop
        brange;
    end
    % extract the dataset to apply the GMM on 
    [test_data, ~] = prep_data(ao_data, brange);
    test_data = unique(test_data);
    % get statistics and class classifications 
    [idx, post_pdf, npost_pdf, d2] = evalPoints(gmm, test_data, 'type', 'mdist');
    mixed_pdf = post_pdf.*npost_pdf;
    gnpost_pdf = post_pdf ./ sum(post_pdf, 1);
    % for each cluster 
    

    for i = 1:nclusters
        
        % get the peak values that correpsond to clsuter i 
        cluster_i_cand = test_data;
        % get the associated distances that correspond to cluster i 
        d2_i_cand = d2(:, i);
        d2_i_cand = d2_i_cand'; % transpose because it looks for max/min along rows 
        prob_i_cand = post_pdf(:, i)';
        mixed_i_cand = mixed_pdf(:, i)';
        % TO DO (POTENTIALLY): 
        %       if we do some normalization (like for probability) remove
        %       probabilities that correspond to identical peak locations 

        % find the idx in the cluster values taht has the best parameters 
        cluster_i_idx_dist = gmm.findCluster(d2_i_cand, 'min');
        cluster_i_idx_prob = gmm.findCluster(prob_i_cand, 'max');
        cluster_i_idx_mixed = gmm.findCluster(mixed_i_cand, 'max');


        % save the value 
        est_rao_dist(brange, i) = cluster_i_cand(cluster_i_idx_dist);
        est_rao_prob(brange, i) = cluster_i_cand(cluster_i_idx_prob);
        est_rao_mixed(brange, i) = cluster_i_cand(cluster_i_idx_mixed);
        
    end 
    
    % save the indices to easily access the destination node info in other
    % matrices 
    thresh = 0.1;
    test_matrix = gnpost_pdf;
    test_matrix = [zeros(size(test_matrix, 1), 1), test_matrix];
    test_matrix(1, 1) = 1;
    value_matrix = d2;
    value_matrix = [zeros(size(value_matrix, 1), 1), value_matrix];
    dijk = Dijkstra();
    dijk = dijk.prep_dijkstra(test_matrix, thresh, value_matrix);
    
    [node_dist, node_routes] = dijk.calc_paths(dijk.start_nodes, dijk.end_nodes);
    
    % need to alter column because first node offsets the matrix by 1
    dijk = dijk.setval('colidx', dijk.colidx - 1);
    pk_vals_routes = test_data(dijk.rowidx(node_routes(:, 2:end)));
       
    % now find the best route using the min distance 
    [best_dist, best_route_idx] = min(node_dist);
    best_node_route = node_routes(best_route_idx, 2:end);
    if length(node_dist) ~= 1
        best_pk_val_route = pk_vals_routes(best_route_idx, :);
    else
        best_pk_val_route = pk_vals_routes;
    end

    est_rao_dijk(brange, :) = best_pk_val_route;
   
    update_matrix = est_rao_dijk;

    if kFpeaks
    
        % kF obs dependent on which peaks you are choosing 
        kFobs = est_rao_dijk(brange, :);

        node_ind = sub2ind(size(d2), dijk.rowidx(best_node_route), ...
            dijk.colidx(best_node_route))';

         % for now just try the d2 values
         R_vec = d2(node_ind);

         all_dist(brange, :) = R_vec;
         kFthresh = 1e6*ones(1, nstates);
         ss_R = diag(min(abs([R_vec; kFthresh])));
    %     %%%% generate Q 
         ss_Q = q*covMatrix;
    %     
         [stateEst(brange, :), kFstatePred, P_pred, res(brange, :)] = ...
             kalmanFilterAlt(ss_A, ss_C, ss_Q, ss_R, ss_S, ...
             kFobs, 'forLoop', kFstatePred, P_pred);

         % for the updating single step mu 
        update_matrix = stateEst;
        
    end
    
    
    if kFmulti
        
        % here we are assuming the matrix is 
        % [p1, p2, p3, p4, u1, u2, u3, u4];
        % kF obs dependent on which peaks you are choosing 
        peakobs = est_rao_dijk(brange, :);
        muobs = reshape(gmm.mu, 1, []);
        kFobs = [peakobs, muobs];
        
        node_ind = sub2ind(size(d2), dijk.rowidx(best_node_route), ...
            dijk.colidx(best_node_route))';
         peakR = d2(node_ind);
         muR = peakR*mu_val; % this might be wrong 
         % for now just try the d2 values
         R_vec = [peakR, muR];

         all_dist(brange, :) = R_vec;
         kFthresh = 1e6*ones(1, nstates);
         ss_R = diag(min(abs([R_vec; kFthresh])));
         ss_Q = q*covMatrix;

         % for the updating single step mu 
         [stateEst(brange, :), kFstatePred, P_pred, res(brange, :)] = ...
             kalmanFilterAlt(ss_A, ss_C, ss_Q, ss_R, ss_S, ...
             kFobs, 'forLoop', kFstatePred, P_pred);update_matrix = est_rao_dijk;
         
         
         stateEst;
         
         % update the mu
         gmm.mu = stateEst(brange, 5:end)';
         est_rao_multi(brange, :) = stateEst(brange, 1:4);
    end
    
    
    % if we pass a certain number of beats 
    if brange > 100 & update & ~kFmulti

            % update mu 
            old_mu = gmm.mu;
            new_pks = update_matrix(brange, :)';
            gmm.mu = old_mu + (new_pks - old_mu)/100;


    end
    
    if brange == bstop

        figure; plot(untesting_beats(:, brange)); hold on; 
        scatter(rao_cand(brange, :), rao_amps(brange, :), 300);
        scatter(est_rao_dist(brange, :), untesting_beats(est_rao_dist(brange, :), brange), 'k')
        scatter(est_rao_dijk(brange, :), untesting_beats(est_rao_dijk(brange, :), brange), 200, '+')
        legend('Beat', 'Peaks', 'Using Just Min on Each Cluster', 'Dijkstra')

    end
    
    all_mu(brange, :) = gmm.mu;
    
    
    %%%% Ok now consider using the kalman filter
    
    
    
    %%%%
end

if ~update && ~kFpeaks && ~kFmulti
    dijknoupdate = est_rao_dijk;
elseif update && ~kFpeaks && ~kFmulti
    dijkupdate = est_rao_dijk;
elseif update && kFpeaks && ~kFmulti
    pkkfupdatekf = stateEst; % kf smoothened pks for updates
    pkkfupdatedi = est_rao_dijk; % dijk pks for updates
elseif update && ~kFpeaks && kFmulti
    multiupdate = est_rao_multi;
end



end
%%
%
if ndim == 1
    lnwid = 2;
    estlnwid = 2;
    f = figure; 
    tiledlayout(6, 1, 'Padding', 'none', 'TileSpacing', 'compact');
    
    
    %
    % dijk no update no kf 
    ax1 = nexttile;
    [~, dijknoupdate_idx] = min(abs(dijknoupdate(20, :) - rao(20)*2));
    dijknoupdate_rmse = rms(rao(1:length(rdict)) - dijknoupdate(:, dijknoupdate_idx)/2);
    plot(dijknoupdate/2);  hold on;
    p1 = plot(dijknoupdate(:, dijknoupdate_idx)/2, 'r' , 'LineWidth', estlnwid);
    p = plot(rao(1:length(rdict)), 'k' , 'LineWidth', lnwid);  
    p.Color(4) = 0.35; ylabel('ms')

    title(['DIJK no \mu update no KF RMSE: ', num2str(dijknoupdate_rmse)])
    
    % dijk update kf plot
    ax2 = nexttile;
    [~, dijkupdate_idx] = min(abs(dijkupdate(20, :) - rao(20)*2));
    dijkupdate_rmse = rms(rao(1:length(rdict)) - dijkupdate(:, dijkupdate_idx)/2);
    plot(dijkupdate/2); hold on;
    p1 = plot(dijkupdate(:, dijkupdate_idx)/2, 'r' , 'LineWidth', estlnwid);
    p = plot(rao(1:length(rdict)), 'k' , 'LineWidth', lnwid); p.Color(4) = 0.35;
    ylabel('ms')

    title(['DIJK \mu update no KF RMSE: ', num2str(dijkupdate_rmse)])
    
    % dijk update with kf kf outputs 
    ax3 = nexttile;

    [~, pkkfupdatekf_idx] = min(abs(pkkfupdatekf(20, :) - rao(20)*2));
    pkkfupdatekf_rmse = rms(rao(1:length(rdict)) - pkkfupdatekf(:, pkkfupdatekf_idx)/2);
    plot(pkkfupdatekf/2); hold on;
    p1 = plot(pkkfupdatekf(:, pkkfupdatekf_idx)/2, 'r' , 'LineWidth', estlnwid);
    p = plot(rao(1:length(rdict)), 'k' , 'LineWidth', lnwid); 
    p.Color(4) = 0.35; ylabel('ms')
    title(['DIJK \mu update using KF RMSE: ', num2str(pkkfupdatekf_rmse)]) 
    
    % dijk update with kf kf outputs and mu adaptations  
    ax4 = nexttile;
    [~, multiupdate_idx] = min(abs(multiupdate(20, :) - rao(20)*2));
    multiupdate_rmse = rms(rao(1:length(rdict)) - multiupdate(:, multiupdate_idx)/2);
    plot(multiupdate/2); hold on;
    p1 = plot(multiupdate(:, multiupdate_idx)/2, 'r' , 'LineWidth', estlnwid);
    p = plot(rao(1:length(rdict)), 'k' , 'LineWidth', lnwid);  p.Color(4) = 0.35;
     ylabel('ms')

    title(['DIJK \mu update using multiKF  RMSE: ', num2str(multiupdate_rmse)]) 
    
    % other
    ax5 = nexttile;
    % calc best idx for sliding template method 
    [~, sltemp_idx] = min(abs(altmethod(20, :) - rao(20)*2));
    sltemp_rmse = rms(rao(1:length(rdict)) - altmethod(:, sltemp_idx)/2);
    plot(altmethod/2); hold on;
    p1 = plot(altmethod(:, sltemp_idx)/2, 'r' , 'LineWidth', estlnwid); 
    p = plot(rao(1:length(rdict)), 'k' , 'LineWidth', lnwid);  p.Color(4) = 0.35;
     ylabel('ms')
    legend([p, p1], 'trueRAO', 'estRAO', 'orientation', 'horizontal')
    title(['Sliding Template RMSE: ', num2str(sltemp_rmse)]) 
    
    ax6 = nexttile;
    %plot(abs(rao(1:length(rdict)) - dijknoupdate(:, dijknoupdate_idx)/2)); hold on;
    %plot(abs(rao(1:length(rdict)) - dijkupdate(:, dijkupdate_idx)))
    plot(abs(rao(1:length(rdict)) - altmethod(:, sltemp_idx)/2)); hold on;
    plot(abs(rao(1:length(rdict)) - pkkfupdatekf(:, pkkfupdatekf_idx)/2)); 
    plot(abs(rao(1:length(rdict)) - multiupdate(:, multiupdate_idx)/2));
    legend('Sliding Template', 'DIJK KF \mu update', 'DIJK multi KF \mu update', 'orientation', 'horizontal')
    title(['Abs error']); xlabel('Beats')  
    linkaxes([ax1, ax2, ax3, ax4, ax5 ax6], 'x');
    sgtitle([num2str(sub), ' ', cln_str])
    savefig(f, [savePath, filesep, num2str(sub), '_', cln_str, '_rao', '.fig'])

    f = figure; 
    %tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'compact');
    ax1 = subplot(2, 1, 1);
    plot(abs(rao(1:length(rdict)) - altmethod(:, sltemp_idx)/2)); hold on;
    plot(abs(rao(1:length(rdict)) - pkkfupdatekf(:, pkkfupdatekf_idx)/2)); hold on;
    plot(abs(rao(1:length(rdict)) - multiupdate(:, multiupdate_idx)/2))
    title(['Abs error']); xlabel('Beats')    
    legend('Sliding Template', 'DIJK KF \mu update', ...
        'DIJK multi KF \mu update')
    ax2 = subplot(2, 1, 2);
    plot(1:length(rdict), rdict);
    ylabel('SNR (dB')
    linkaxes([ax1, ax2], 'x');
    savefig(f, [savePath, filesep, num2str(sub), '_', cln_str, '_error', '.fig'])
    
    stats = struct; 
    stats.sub = sub;
    stats.dijknoupdate = dijknoupdate;
    stats.dijkupdate = dijkupdate;
    stats.pkkfupdatekf = pkkfupdatekf;
    stats.pkkfupdatedi = pkkfupdatedi;
    stats.multiupdate = multiupdate;
    stats.slidetemp = altmethod; 
    stats.rao = rao(1:length(rdict)); 
    stats.rdict = rdict;
    stats.rmse = [dijknoupdate_rmse, dijkupdate_rmse, pkkfupdatekf_rmse, ...
         multiupdate_rmse, sltemp_rmse];
    save([savePath, filesep, num2str(sub), '_', cln_str, '.mat'], 'stats')

    %linkaxes([ax1, ax2, ax3, ax4, ax5 ], 'x');
    
    
end


end
close all
end
%%
sub = 1;
load([savePath, filesep, num2str(sub), '_', cln_str, '.mat'])
sublist = 1:6;
% find out how many snr levels there are 
snrlist = (unique(stats.rdict));

% snr  levels [levels, all, clean]
snr_rmse_list = zeros(length(snrlist)+2, length(sublist), 5);


for sub = sublist
    load([savePath, filesep, num2str(sub), '_', 'cln', '.mat'])
    rao = stats.rao;
    noupnokf = findRAO(stats.dijknoupdate, stats.rao)/2;
    upnokf = findRAO(stats.dijkupdate, stats.rao)/2;
    upkf = findRAO(stats.pkkfupdatekf, stats.rao)/2;
    upmultikf = findRAO(stats.multiupdate, stats.rao)/2;
    sldtemp = findRAO(stats.slidetemp, stats.rao)/2;
    
    load([savePath, filesep, num2str(sub), '_', 'noisy', '.mat'])
    noupnokf_noisy = findRAO(stats.dijknoupdate, stats.rao)/2;
    upnokf_noisy = findRAO(stats.dijkupdate, stats.rao)/2;
    upkf_noisy = findRAO(stats.pkkfupdatekf, stats.rao)/2;
    upmultikf_noisy = findRAO(stats.multiupdate, stats.rao)/2;
    sldtemp_noisy = findRAO(stats.slidetemp, stats.rao)/2;
    rdict = stats.rdict;

    snr_rmse_list(:, sub, 1) = [rmseSeg(rao, noupnokf_noisy, rdict, -20); ...
        rmseSeg(rao, noupnokf_noisy, rdict, -15); rmseSeg(rao, noupnokf_noisy, rdict, -10); ...
        rmseSeg(rao, noupnokf_noisy, rdict, -5); rmseSeg(rao, noupnokf_noisy, rdict, 0); ...
        rmseSeg(rao, noupnokf_noisy, rdict, 5); rmseSeg(rao, noupnokf_noisy, rdict, 10); ...
        rmseSeg(rao, noupnokf_noisy, rdict, 'all'); rmseSeg(rao, noupnokf, rdict, 'all')];
    

    snr_rmse_list(:, sub, 2) = [rmseSeg(rao, upnokf_noisy, rdict, -20); ...
        rmseSeg(rao, upnokf_noisy, rdict, -15); rmseSeg(rao, upnokf_noisy, rdict, -10); ...
        rmseSeg(rao, upnokf_noisy, rdict, -5); rmseSeg(rao, upnokf_noisy, rdict, 0); ...
        rmseSeg(rao, upnokf_noisy, rdict, 5); rmseSeg(rao, upnokf_noisy, rdict, 10); ...
        rmseSeg(rao, upnokf_noisy, rdict, 'all'); rmseSeg(rao, upnokf, rdict, 'all')];
    

    snr_rmse_list(:, sub, 3) = [rmseSeg(rao, upkf_noisy, rdict, -20); ...
        rmseSeg(rao, upkf_noisy, rdict, -15); rmseSeg(rao, upkf_noisy, rdict, -10); ...
        rmseSeg(rao, upkf_noisy, rdict, -5); rmseSeg(rao, upkf_noisy, rdict, 0); ...
        rmseSeg(rao, upkf_noisy, rdict, 5); rmseSeg(rao, upkf_noisy, rdict, 10); ...
        rmseSeg(rao, upkf_noisy, rdict, 'all'); rmseSeg(rao, upkf, rdict, 'all')];
    
    snr_rmse_list(:, sub, 4) = [rmseSeg(rao, upmultikf_noisy, rdict, -20); ...
        rmseSeg(rao, upmultikf_noisy, rdict, -15); rmseSeg(rao, upmultikf_noisy, rdict, -10); ...
        rmseSeg(rao, upmultikf_noisy, rdict, -5); rmseSeg(rao, upmultikf_noisy, rdict, 0); ...
        rmseSeg(rao, upmultikf_noisy, rdict, 5); rmseSeg(rao, upmultikf_noisy, rdict, 10); ...
        rmseSeg(rao, upmultikf_noisy, rdict, 'all'); rmseSeg(rao, upmultikf, rdict, 'all')];
    
    snr_rmse_list(:, sub, 5) = [rmseSeg(rao, sldtemp_noisy, rdict, -20); ...
        rmseSeg(rao, sldtemp_noisy, rdict, -15); rmseSeg(rao, sldtemp_noisy, rdict, -10); ...
        rmseSeg(rao, sldtemp_noisy, rdict, -5); rmseSeg(rao, sldtemp_noisy, rdict, 0); ...
        rmseSeg(rao, sldtemp_noisy, rdict, 5); rmseSeg(rao, sldtemp_noisy, rdict, 10); ...
        rmseSeg(rao, sldtemp_noisy, rdict, 'all'); rmseSeg(rao, sldtemp, rdict, 'all')];

end


snr_rmse_list = permute(snr_rmse_list, [2, 1, 3]);

subcolors = [0, 0.4470, 0.7410;
    0.85, 0.325, 0.098;
    0.929, 0.694, 0.125;
    0.494, 0.184, 0.556;
    0.466, 0.674, 0.188;
    0.301, 0.745, 0.933];
[nsubs, nsnrs, nmethods] = size(snr_rmse_list);

% number of box plots 
nboxplots = nsnrs * nmethods;
snr_labels = {'-20', '-15', '-10', '-5', '0', '5', '10', 'N', 'C'};
method_labels =  {'nUnK', 'UnK', 'UK', 'UMK', 'S'};
% what elements you want i.e what should be each of the columns of the
% matrix 
type = 'method';
if strcmp('snr', type)
    % number of groups (assume groups are snr: 9 or groups are methods: 5
    nelems = nsnrs;  ngroups = nmethods; rmse_list = permute(snr_rmse_list, [1, 3, 2]);
    primary_labels = snr_labels;
    group_labels = method_labels;
elseif strcmp(type, 'method')
    nelems = nmethods; ngroups = nsnrs; rmse_list = (snr_rmse_list);
    primary_labels = method_labels;
    group_labels = snr_labels;
end
gspace = 1;
data_cell = cell(1, nelems);
for i = 1:size(data_cell, 2)
   [row, col] = ind2sub([size(rmse_list, 1), size(rmse_list, 2)], find(squeeze(rmse_list(:, :, i)) > inf));
   rmse_list(row, col, i) = nan;
   data_cell{i} = rmse_list(:, :, i); 
   data_cell{i}(data_cell{i} > 20) = nan;
end


pos_vec = [];
alt_rmse_list = [];
for i = 1:ngroups
    pos_vec =  [pos_vec, (1:nelems) + (i-1)*(nelems+gspace)];
    
end
pos_vec = repmat(pos_vec, 6, 1);
pos_vec = reshape(pos_vec, 1, []);

for i = 1:ngroups
    alt_rmse_list = [alt_rmse_list; reshape(rmse_list(:, i, :), [], 1)];
end

% alt_rmse_list = reshape(rmse_list, nsubs, nelems*ngroups);
% alt_rmse_list = reshape(alt_rmse_list, [], 1)
colors = repmat(subcolors', 1, ngroups*nelems);
f = figure; 
scatter(pos_vec, alt_rmse_list, 50, colors', 'filled'); hold on;
s = boxplotGroup(data_cell, 'PrimaryLabels', primary_labels, 'SecondaryLabels', group_labels, ...
    'InterGroupSpace', gspace);
title('RMSE Comparisons (Peak Switching not Considered)'); 
ylabel('RMSE (ms)')
savefig(f, [savePath, filesep, 'aggregate_rmse.fig'])

%%
lnwid = 2;
estlnwid = 2;
f = figure; 
tiledlayout(3, 1);

load([savePath, filesep, num2str(2), '_', 'noisy', '.mat'])
noupnokf_noisy = findRAO(stats.dijknoupdate, stats.rao)/2;
upnokf_noisy = findRAO(stats.dijkupdate, stats.rao)/2;
upkf_noisy = findRAO(stats.pkkfupdatekf, stats.rao)/2;
upmultikf_noisy = findRAO(stats.multiupdate, stats.rao)/2;
sldtemp_noisy = findRAO(stats.slidetemp, stats.rao)/2;
rao = stats.rao;
rdict = stats.rdict
%
% dijk no update no kf 
ax1 = nexttile;
plot(stats.slidetemp(:, 1)/2); hold on;
p1 = plot(sldtemp_noisy, 'r', 'LineWidth', estlnwid); 
p = plot(rao, 'k' , 'LineWidth', lnwid); p.Color(4) =1;
title(['Sliding Template RMSE: ', num2str(rmseSeg(rao, sldtemp_noisy, rdict, 'all'))])

ax2 = nexttile; 
plot(stats.dijkupdate(:, 1)/2); hold on;
p1 = plot(upnokf_noisy, 'r', 'LineWidth', estlnwid); hold on;
p = plot(rao, 'k' , 'LineWidth', lnwid); p.Color(4) =1;
title(['Dijkstra single step \mu update no KF RMSE: ', num2str(rmseSeg(rao, upnokf_noisy, rdict, 'all'))])

a3 = nexttile;
p2 = plot(stats.multiupdate(:, 1)/2); hold on;
p1 = plot(upmultikf_noisy, 'r', 'LineWidth', estlnwid); hold on;
p = plot(rao, 'k' , 'LineWidth', lnwid); p.Color(4) =1;
legend([p, p1, p2], {'True AO', 'Est AO', 'First Peak'})
title(['Dijkstra Multi Kalman \mu update RMSE: ', num2str(rmseSeg(rao, upmultikf_noisy, rdict, 'all'))])

ylabel('R-Peak Distance (ms)');
xlabel('Beats')
  
%% functions 


function noise = prod_noise(sig, type, seed)
% 
rng(seed)
[sigX, sigY] = size(sig);
% take the square root?
if strcmp(type, 'gauss')
    noise = randn(sigX, sigY);
elseif strcmp(type, 'chirp')
    fs = 2000;
else 
    error('input gauss')
end
end


function estrao = findRAO(candidates, truerao)
    % candidates Nx nclusters 
    [~, estrao_idx] = min(abs(candidates(20, :) - truerao(20)*2));
    estrao = candidates(:, estrao_idx);
end

function rmse_seg = rmseSeg(estrao, truerao, rdict, delineator)
    if isstr(delineator)
        if strcmp(delineator, 'all')
           rmse_seg = rms(truerao - estrao); 
        else
            error('delineator should be all')
        end
    
    else
        rmse_seg = rms(truerao(rdict == delineator) - estrao(rdict == delineator));
    end

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
function rao_cand = findAO(sepbeats, num_cand, fs)
    %%%%%%%%%%%%%% Calculate top 10 peaks in each beat within 250ms %%%%%%%%%%
    rao_cand = zeros(size(sepbeats, 2), num_cand);

    for beat = 1:size(sepbeats, 2)


        temp_beat = sepbeats(:, beat);
        [pks, locs] = findpeaks(temp_beat(1:250*fs/1000), 'MinPeakDistance', 10*fs/1000);
        if length(locs) < num_cand
            locs = [locs; ones(num_cand - length(locs), 1)*locs(end)];
            pks = [pks; ones(num_cand - length(pks), 1)*pks(end)];
        end
        [~, idx] = sort(pks, 'descend');

        rao_cand(beat, :) = sort(locs(idx(1:num_cand)), 'ascend');


    end
end

function [datavec, groupvec, org_dim] = prep_data(data, range)
% here we assume data is a N x ncluster x ndim dataset 
% makes a group vector Mxndim where M <= N 
% makes a data vector Mxndim  where M <= N
% range specifies the region to concatenate300

data = data(range, :, :);
[N, nclusters, ndims] = size(data);

group = repmat(1:nclusters, N, 1);

datavec = zeros(N*nclusters, ndims);
groupvec = reshape(group, [], 1);

for i = 1:ndims
    datavec(:, i) = reshape(data(:, :, i), [], 1);
end

org_dim = [N, nclusters, ndims];
end
function [data] = revert_data(data, org_dim);

reshape(data, org_dim);

end

function gmmMovie(signals, models, varargin)

% -------------------------------------------------------------------------
% This function creates an animation of a heartbeat-separated signal with
% optional feature points overlaid.
%
% Arguments (required)
% - signals     [MxN]   N signal vectors of length M for visualization
%
% Arguments (optional)
% - 'features'  [NxF]   Feature points for visualization (F per heartbeat)
% - 'start'             Sample at which to start animation
% - 'stop'              Sample at which to stop animation
% - 'speed'             Animation speed ('slow', 'med', or 'fast')
% - 'manual'            Change beat examined using arrow keys or mouse clicks
% 
% Usage: heartbeatMovie(signals, varargin)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'features'); features = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'start'); start = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'stop'); stop = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'speed'); speed = varargin{arg + 1};
        elseif strcmp(varargin{arg}, 'manual'); manual = true;
        end
    end
end

% Get the number of signals and features
numSignals = length(signals);

% Set defaults for optional arguments
if ~exist('start', 'var'); start = 1; end
if ~exist('stop', 'var'); stop = numSignals; end
if ~exist('speed', 'var'); speed = 'med'; end
if ~exist('manual', 'var'); manual = false; end

% Set the speed parameter tau
switch speed
    case 'slow'
        tau = 1;
    case 'med'
        tau = 0.1;
    case 'fast'
        tau = 0.01;
    otherwise
        tau = 0.1;
end
[N, ndim] = size(signals{1});

% figure for plotting signals 
figure


    for i = start:stop
        pointsz = 20;
        data = signals{i};
        gmm = models{i};
        if ndim == 2
            % Plot the current signal
            scatter(data(:,1),data(:,2),pointsz,'.') 
            hold on
            % https://www.mathworks.com/help/stats/gmdistribution.pdf.html
            gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gmm.model,[x0 y0]),x,y);
            plotdim = [min(data(:, 1)), max(data(:, 1)), min(data(:, 2)), max(data(:, 2))];
            fcontour(gmPDF,plotdim) 
        elseif ndim == 1
            gmPDF = @(x) arrayfun(@(x0) pdf(gmm.model,[x0]), x);
            scatter(data, gmPDF(data)); hold on;
            plotdim = [min(data), max(data)];
            fplot(gmPDF,plotdim)   
        end

        % Display figure title
        title("Segment " + string(i) + " of " + string(stop - start))


        % Pause and reset figure
        pause(tau); clf('reset');


    end

end





function handles = boxplotGroup(varargin)
% BOXPLOTGROUP groups boxplots together with horizontal space between groups.
%   boxplotGroup(x) receives a 1xm cell array where each element is a matrix with
%   n columns and produces n groups of boxplot boxes with m boxes per group.
%
%   boxplotGroup(ax,x,___) specifies the axis handle, otherwise current axis is used.
%
%   boxplotGroup(___,'interGroupSpace',d) separates groups by d units along the x axis
%   where d is a positive, scalar integer (default = 1)
%
%   boxplotGroup(___,'groupLines', true) adds vertical divider lines between groups
%   (requires >=r2018b).
%
%   boxplotGroup(___,'primaryLabels', c) specifies the x tick label for each boxplot.
%   c is a string array or cell array of characters and must have one element per
%   box or one element per group-member. When undefined or when c is an empty cell {},
%   the x-axis is labeled with default x-tick labels.
%
%   boxplotGroup(___,'secondaryLabels', s) specifies the group labels for the boxplot
%   groups.  s is a string array or cell array of characters and must have one element
%   per group (see 'groupLabelType'). Ignored when s is an empty cell {}.
%
%   boxplotGroup(___,'groupLabelType', str) specifies how to label the groups by one of
%   the following options.
%    * 'horizontal': Group labels will be centered under the primary labels using a 2nd
%       invisible axis underlying the main axis (not supported in uifigures). To remove
%       the primary labels and only show secondary labels, set primary labels to empty
%       cell-strings (e.g. {'','',''}) or strings without characters (e.g. ["" "" ""]).
%    * 'vertical': Group labels will be vertical, between groups (requires Matlab >=2018b)
%    * 'both': Both methods will be used.
%
%   boxplotGroup(___, 'PARAM1', val1, 'PARAM2, val2, ...) sends optional name/value pairs
%   to the boxplot() function. Accepted parameters are BoxStyle, Colors, MedianStyle,
%   Notch, OutlierSize, PlotStyle, Symbol, Widths, DataLim, ExtremeMode, Jitter, and Whisker.
%   See boxplot documentation for details.
%
%   boxplotGroup(___, 'Colors', ___, 'GroupType', type) determines how to apply
%   colors to the groups.  'Colors' is a property of boxplots (see boxplot documentation).
%   When the colors value specifies multiple colors, the 'GroupType' determines how
%   the colors are distributed based on the following two options.
%    * 'betweenGroups' assigns color n to the n^th boxplot within each group (default).
%    * 'withinGroups' assigns color n to all boxplots within the n^th group.
%
%   h = boxplotGroup(___) outputs a structure of graphics handles.
%
% NOTE: If you're working with a grouping variable 'g', use the syntax boxplot(x,g) along
%   with the "Group Appearance" options described in Matlab's boxplot() documentation.
%   https://www.mathworks.com/help/stats/boxplot.html#d118e146984
%
% EXAMPLES:
% data = {rand(100,4), rand(20,4)*.8, rand(1000,4)*1.2};
%
% Required inputs
%   boxplotGroup(data)
%
% Set space between groups
%   boxplotGroup(data, 'interGroupSpace', 3)
%
% Specify labels and draw divider line
%   boxplotGroup(data, 'groupLines', true, 'PrimaryLabels', {'a' 'b' 'c'},...
%       'SecondaryLabels', {'Lancaster', 'Cincinnati', 'Sofia', 'Rochester'})
%
% Label groups with vertical lables
%   boxplotGroup(data, 'PrimaryLabels', {'a' 'b' 'c'}, 'SecondaryLabels', ...
%       {'Lancaster', 'Cincinnati', 'Sofia', 'Rochester'}, 'groupLabelType', 'vertical')
%
% Pass additional boxplot properties
%   boxplotGroup(data, 'PrimaryLabels', {'a' 'b' 'c'}, 'SecondaryLabels', ...
%       {'Lancaster', 'Cincinnati', 'Sofia', 'Rochester'}, 'groupLabelType', 'vertical', ...
%       'BoxStyle', 'filled', 'PlotStyle', 'Compact')
%
%
% Contact adam.danz@gmail.com for questions, bugs, suggestions, and high-fives.
% Copyright (c) 2020, Adam Danz  adam.danz@gmail.com
% All rights reserved
% Source: https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup

% Changes history
% 200306 - v1.0.0 first release.
% 200308 - v1.1.0 Added recommendation to use boxplot() with grouping variable.
%                 Added axis handle as input to boxplot() call. Linkaxes changed
%                 from x to xy. Set axis2.Units to axis.Units.  Using linkprop
%                 to link position etc of main axis and axis2. Added DeleteFcn
%                 to main axis. Disabled toolbar for axis2. Added listener to
%                 resize axis2 when main axis is resized. Changes to help section.
% 200309 - v1.2.0 When 2nd axis is added, main axis is set to current axis.
% 200309 - v1.2.1 Suppress linkprops() and changes to toolbar suppression to work
%                 with versions prior to r2018b.
% 200309 - v1.2.2 Instead of creating new axis, default axis is gca().
% 210427 - v2.0.0 oncleanup returns hold state instead of conditional.  Added GroupType
%                 option and colorexpansion. Suppresses output unless requested.
%                 Checks matlab vs with xline(). Removed listener, storing hlink in axis.
%                 boxplot name-val arg check. Removing boxplot placeholders. XTicks now auto
%                 if labels aren't provided. Outputs now include boxplotGroup; vertical
%                 labels now the same fontsize and weight as axis font; Primary and secondary
%                 labels can be empty cell to ignore. Secondary labels now match ax1 font size,
%                 weight and name.

%% Check for axis handle in first input
if ~isempty(varargin) && ~isempty(varargin{1}) && isobject(varargin{1}(1)) % [3]
    if isgraphics(varargin{1}(1), 'axes')
        % first input is an axis
        h.axis = varargin{1} ;
        varargin(1) = [];
    else
        error('MATLAB:hg:InvalidHandle', 'Invalid handle')
    end
else
    h.axis = [];
end

%% Parse inputs
p = inputParser();
p.FunctionName = mfilename;
p.KeepUnmatched = true;	%accept additional parameter value inputs (passed to boxplot())
addRequired(p, 'x', @(x)validateattributes(x,{'cell'},{'row','nonempty'}))
addParameter(p, 'interGroupSpace', 1, @(x)validateattributes(x,{'double'},{'scalar','integer'}))
addParameter(p, 'primarylabels', [], @(x)validateattributes(x,{'string','cell'},{}))
addParameter(p, 'secondarylabels', [], @(x)validateattributes(x,{'string','cell'},{}))
addParameter(p, 'groupLines', false, @(x)validateattributes(x,{'logical','double'},{'binary'}))
addParameter(p, 'groupLabelType', 'Horizontal', @(x)ischar(validatestring(lower(x),{'vertical','horizontal','both'})))
addParameter(p, 'GroupType', 'betweenGroups', @(x)ischar(validatestring(lower(x),{'betweengroups','withingroups'})))
parse(p,varargin{:})

% Prepare the unmatched boxplot() parameters.
% If a param is passed that isn't accepted by boxplot(), an error is thrown from boxplot() function.
unmatchNameVal = reshape([fieldnames(p.Unmatched)'; struct2cell(p.Unmatched)'], 1, []);

% Check boxplot name-value parameters; group params, Position, and labels are not accepted.
supportedParams = {'BoxStyle','Colors','MedianStyle','Notch','OutlierSize','PlotStyle','Symbol','Widths', ...
    'DataLim','ExtremeMode','Jitter','Whisker'};
argOK = arrayfun(@(i)any(strncmpi(unmatchNameVal{i},supportedParams,numel(unmatchNameVal{i}))),...
    1:2:numel(unmatchNameVal)); % look for partial match
assert(all(argOK),'Parameter(s) not accepted in %s: [%s].', ...
    mfilename, strjoin(unmatchNameVal(find(~argOK)*2-1),', '))

% Check that each element of x is a matrix
assert(all(cellfun(@ismatrix, p.Results.x)), 'All elements of the cell array ''x'' must be a matrix.')
% Check that each matrix contains the same number of columns.
assert(numel(unique(cellfun(@(m)size(m,2),p.Results.x))) == 1, ...
    ['All elements of the cell array ''x'' must contain the same number of columns. '...
    'Pad the matricies that contain fewer columns with NaN values.']);

nargoutchk(0,1)

%% Compute horizontal spacing & check labels
nGroups = size(p.Results.x{1},2);       % number of columns of data / number of groups
nMembers = numel(p.Results.x);          % number of members per group
maxX = ((nMembers + p.Results.interGroupSpace) * nGroups) - p.Results.interGroupSpace;
xInterval = nMembers + p.Results.interGroupSpace;

% Check that labels (if any) are the right size
% PrimaryLabels: either 1 per group-member or 1 for each bar
if ~isempty(p.Results.primarylabels)
    assert(ismember(numel(p.Results.primarylabels),[nMembers, nMembers*nGroups]), ...
        sprintf(['The number of primary labels must equal either the number of bars per group (%d) '...
        'or the number of total bars (%d).'], nMembers, nMembers*nGroups))
end
% SecondaryLabels: 1 per group
if ~isempty(p.Results.secondarylabels)
    assert(isequal(numel(p.Results.secondarylabels),nGroups), ...
        sprintf('The number of secondary labels must equal either the number groups (%d).',nGroups))
end

% If all primary labels are empty chars do not add the newline to secondary labels.
if ~isempty(p.Results.primarylabels) &&  all(cellfun(@isempty,cellstr(p.Results.primarylabels)))
    horizSecondaryLabelAddon = '';
else
    horizSecondaryLabelAddon = '\newline';
end

%% Set colors
% Assumes ColorGroup property is not specified.
colorsIdx = strcmpi('Colors',unmatchNameVal);
if any(colorsIdx)
    cvalIdx = find(colorsIdx,1,'first')+1;
    if isempty(unmatchNameVal{cvalIdx})
        % Colors val is empty; remove Colors name-val pair
        unmatchNameVal(cvalIdx-[1,0]) = [];
    else
        unmatchNameVal{cvalIdx} = colorexpansion(unmatchNameVal{cvalIdx}, p, nGroups, nMembers);
    end
end

%% Do plotting
if isempty(h.axis)
    h.axis = gca();
end
h.figure = ancestor(h.axis,'figure');
isTiledLayout = strcmpi(h.axis.Parent.Type,'tiledlayout');
if isTiledLayout % [6]
    origTLOState = warning('query', 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout');
    TLOcleanup = onCleanup(@()warning(origTLOState));
    warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')
end

% Store pre-existing boxplot object handles
bptag = 'boxplot'; % tag Matlab assigns to bp group
bpobjPre = findobj(h.axis,'tag',bptag);

originalHoldStatus = ishold(h.axis);
holdStates = {'off','on'};
returnHoldState = onCleanup(@()hold(h.axis,holdStates{originalHoldStatus+1}));
hold(h.axis, 'on')

x = cell(1,nMembers);
existingTextObjs = findobj(h.axis,'Type','Text');
for i = 1:nMembers
    x{i} = i : xInterval : maxX;
    temp = nan(size(p.Results.x{i},1), max(x{i}));
    temp(:,x{i}) = p.Results.x{i};
    boxplot(h.axis, temp, unmatchNameVal{:})
end

% Remove dummy boxplots placeholders
bpobjNew = findobj(h.axis,'tag',bptag);
bpobjNew(ismember(bpobjNew, bpobjPre)) = [];
for g = 1:numel(bpobjNew)
    tags = unique(get(bpobjNew(g).Children,'tag'),'stable');
    tags(cellfun(@isempty,tags)) = [];
    for j = 1:numel(tags)
        obj = findobj(bpobjNew(g),'tag',tags{j});
        obj(~isprop(obj,'YData')) = [];
        YData = get(obj,'YData');
        if ~iscell(YData)
            YData = {YData};
        end
        isDummy = cellfun(@(c)all(isnan(c),2),YData);
        delete(obj(isDummy))
    end
end

axis(h.axis, 'tight')
limGap = (p.Results.interGroupSpace+1)/2;
set(h.axis,'XTickMode','Auto','XTickLabelMode','Auto','xlim',[1-limGap, maxX+limGap]) %[1]
yl = ylim(h.axis);
ylim(h.axis, yl + [-range(yl)*.05, range(yl)*.05])

% Remove boxplot's text-tics [1]
allTextObjs = findobj(h.axis,'Type','Text');
isBoxplotText = ~ismember(allTextObjs,existingTextObjs);
set(allTextObjs(isBoxplotText), 'String','','Visible','off')

% Set primary labels if provided
if ~isempty(p.Results.primarylabels)
    h.axis.XTick = sort([x{:}]);
    h.axis.XTickLabel = p.Results.primarylabels;
end
% Set secondary labels if provided
vertLinesDrawn = false;
groupLabelType = p.Results.groupLabelType;
if ~isempty(p.Results.secondarylabels)
    if any(strcmpi(groupLabelType, {'horizontal','both'}))
        % Try to detect figure type [4]
        if verLessThan('Matlab','9.0')      %version < 16a (release of uifigs)
            isuifig = @(~)false;
        elseif verLessThan('Matlab','9.5')  % 16a <= version < 18b
            isuifig = @(h)~isempty(matlab.ui.internal.dialog.DialogHelper.getFigureID(h));
        else                                % version >= 18b (written in r21a)
            isuifig = @(h)matlab.ui.internal.isUIFigure(h) && ~isprop(h,'LiveEditorRunTimeFigure');
        end
        isUIFigure = isuifig(h.figure);
        if isUIFigure
            groupLabelType = 'vertical';
            warning('BOXPLOTGRP:uifig','''Horizontal'' GroupLabelType is not supported with UIFIgures. GroupLabelType was changed to ''vertical''.')
        else
            % Tick label rotation must be 0 if using both primary & secondary horizontal labels
            h.axis.XAxis.TickLabelRotation = 0;
            % Compute x position of secondary labels
            if isa(h.axis,'matlab.ui.control.UIAxes')
                axFcn = @uiaxes;
            else
                axFcn = @axes;
            end
            if verLessThan('Matlab','9.8') %r2020a
                posProp = 'Position';
            else
                posProp = 'InnerPosition';
            end
            secondaryX = (nMembers : nMembers + p.Results.interGroupSpace : maxX) - (nMembers-1)/2;
            secondaryLabels = strcat(horizSecondaryLabelAddon,p.Results.secondarylabels); %[2]
            h.axis2 = axFcn(h.figure,'Units',h.axis.Units, 'OuterPosition', h.axis.OuterPosition, ...
                'ActivePositionProperty', h.axis.ActivePositionProperty,'xlim', h.axis.XLim, ...
                'TickLength', [0 0], 'ytick', [], 'Color', 'none', 'XTick', secondaryX, ...
                'TickLabelInterpreter','tex','XTickLabel', secondaryLabels,'HitTest','off',...
                'XTickLabelRotation',0,'box','off','FontSize',h.axis.FontSize,...
                'FontWeight',h.axis.FontWeight,'FontName',h.axis.FontName);
            h.axis.(posProp)([2,4]) = h.axis2.(posProp)([2,4]); % make room in original axes for 2ndary labels.
            h.axis2.(posProp)([1,3]) = h.axis.(posProp)([1,3]); % let original axis control lateral placement
            h.axis2.UserData.hlink = linkprop([h.axis, h.axis2],...
                {'Units',posProp,'ActivePositionProperty','Parent'}); % [5]
            linkaxes([h.axis, h.axis2], 'xy')
            if ~isUIFigure % [4]
                uistack(h.axis2, 'down')
            end
            if isprop(h.axis2, 'Toolbar')
                h.axis2.Toolbar.Visible = 'off'; % ver >= r2018b
            end
            h.axis2.XRuler.Axle.Visible = 'off';
            h.axis2.YRuler.Axle.Visible = 'off';
            h.axis.DeleteFcn = @(~,~)delete(h.axis2); % Delete axis2 if main axis is deleted
            set(h.figure,'CurrentAxes',h.axis)
        end
    end
    if any(strcmpi(groupLabelType, {'vertical','both'})) && ~verLessThan('Matlab','9.5') % r18b
        spaces = setdiff(1-p.Results.interGroupSpace : maxX, [x{:}]);
        endSpaceIdx = [diff(spaces),2] > 1;
        midSpace = spaces(endSpaceIdx) - (p.Results.interGroupSpace-1)/2;
        h.xline = arrayfun(@(x)xline(h.axis, x,'FontSize',h.axis.FontSize,...
            'FontWeight',h.axis.FontWeight,'FontName',h.axis.FontName),midSpace);
        set(h.xline(:), {'Label'}, cellstr(p.Results.secondarylabels(:))) % cellstr in case lbls are str
        vertLinesDrawn = true;
    end
end

% Draw vertical lines if requested and if they don't already exist.
if p.Results.groupLines && ~vertLinesDrawn && ~verLessThan('Matlab','9.5') %r18b
    spaces = setdiff(1:maxX+p.Results.interGroupSpace, [x{:}]);
    endSpaceIdx = [diff(spaces),2] > 1;
    midSpace = spaces(endSpaceIdx) - (p.Results.interGroupSpace-1)/2;
    h.xline = arrayfun(@(x)xline(h.axis, x,'-k'),midSpace);
end

clear('returnHoldState','TLOcleanup')

%% Return output only if requested
if nargout>0
    % Get and organize new boxplot groups
    bpobjPost = findobj(h.axis,'tag',bptag);
    h.boxplotGroup = bpobjPost(~ismember(bpobjPost, bpobjPre));
    handles = h;
end

function c = colorexpansion(colors, p, nGroups, nMembers)
% colors is color data. As of r2021a, boxplot 'Colors' can be RGB triplet/matrix
%   char vec, or string scalar of chars ("rgb"). Long color names is not accepted
%   by boxplot. 'colors' cannot be empty for this function.
% c: if 'colors' specifies more than 1 color, c is the color scheme expanded according
%   to GroupType. Otherwise, c is the same as colors.
% Other inputs defined in main func.
if isnumeric(colors) && size(colors,1)>1
    basecolors = colors;
    
elseif (ischar(colors) || isa(colors,'string')) && numel(char(colors))>1
    basecolors = char(colors);
    basecolors = basecolors(:); % col vec
    
else
    % If colors is not numeric, char, or string let boxplot throw the error.
    % If colors specifies only 1 color, copy colors to output.
    c = colors;
    return
end
isBetweenGroups = strcmpi('betweenGroups', p.Results.GroupType);
n = size(basecolors,1);
getRowIdx = @(n,m)round(mod(1:n,m+1E-08));
if isBetweenGroups
    % The first nMembers of colors will be used
    % Let boxplot do the expansion.
    rowNum = getRowIdx(nMembers,n);
    c = [basecolors(rowNum,:);repmat(basecolors(1,:),p.Results.interGroupSpace,1)];
else
    % The first nGroups colors will be used
    rowNum = getRowIdx(nGroups,n);
    c = repelem(basecolors(rowNum,:),nMembers+p.Results.interGroupSpace,1);
end
if ischar(c)
    c = c';
end


end
end



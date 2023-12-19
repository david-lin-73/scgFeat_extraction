% This script loads in test dat and shows how the sortPeaks algorithm
% is applied to get ao and ac points. It is assumed that the scg signal is
% already filtered and segmented into beats using the ECG. scg_beats should
% be a matrix that is samples/beat x num beats in size
% 
% Recommended parameters to adjust are: 'axrange', 'nclusters', 'influence', 
% 'promconst', 'brange', and findextrema. 
% The other parameters can more or less be set to their default values

% load data
data = load('scg_beats.mat');
scg_beats = data.scg_beats; % extract scg beat matrix [samples/beat x num_beats]
fs = data.fs;               % sampling frequency
findextrema = true;         % automatically find peaks and valleys to track
verbose = true;  

ao_range = [1, 250];        % ao range (in ms)
ac_range = [250, 500];      % ac range (in ms)

disp("Finding AO")
ao = sortPeaks(scg_beats, fs, 'axrange', ao_range, 'verbose', verbose, 'findextrema', findextrema);

disp("Finding AC")
ac = sortPeaks(scg_beats, fs, 'axrange', ac_range, ...
    'nclusters', 3, 'verbose', verbose, 'findextrema', findextrema);


% get most reasonable ao/ac based on standard deviation
[~, best_ao_idx] = min(std(ao));
[~, best_ac_idx] = min(std(ac));
best_ao = ao(:, best_ao_idx);
best_ac = ac(:, best_ac_idx);

% set colors 
cand_color = [0 0.4470 0.7410];
ac_color = [0.4660 0.6740 0.1880];
ao_color = [0.8500 0.3250 0.0980];

% plot heatmap with ao/ac candidates and best ao/ac shown
f = figure; 
f.Position = [675 453 1500 750];
subplot(3, 1, [1, 2])
imagesc(normalize(scg_beats)); c= contrast(normalize(scg_beats)); colormap(c);
hold on; p1 = plot(ao, 'Color', cand_color);  plot(ac, 'Color', cand_color);
p2 = plot(best_ao, 'Color', ao_color); p3 = plot(best_ac, 'Color', ac_color);
title('SCG Beat Heatmap')
xlabel('Beats'); ylabel('Time from R-peak [samps]')
legend([p1(1), p2, p3], {'AO/AC Candidates', 'AO', 'AC'})

% plot scg beats with ao/ac candidates and best ao/ac shown
beats = 1:100;                                          % beats to plot
ao_amps = general.indexMatrix(scg_beats, round(ao));    % ao amplitudes
ac_amps = general.indexMatrix(scg_beats, round(ac));    % ac amplitudes
subplot(3, 1, 3)
p = plot(scg_beats(:, beats), 'k'); hold on;
scatter(ao(beats, :), ao_amps(beats, :), [], cand_color);
s1 = scatter(ac(beats, :), ac_amps(beats, :), [], cand_color);
s2 = scatter(ao(beats, best_ao_idx), ao_amps(beats, best_ao_idx), [], ao_color);
s3 = scatter(ac(beats, best_ac_idx), ac_amps(beats, best_ac_idx), [], ac_color);
title("SCG Beats " + num2str(beats(1)) + " to " +  num2str(beats(end)))
xlabel('Time from R-peak [samps]'); ylabel('SCG Amplitude')
legend([s1(1), s2, s3], {'AO/AC Candidates', 'AO', 'AC'})


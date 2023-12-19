% Ok the goal here is to just give a set of SCG beats and run the algorithm
% and show the plots 


data = load('scg_beats.mat');
scg_beats = data.scg_beats; % extract scg beat matrix [samples/beat x num_beats]
fs = data.fs;               % sampling frequency


% Try runnign the algorithm on the stress data 
path = "/home/dlin/Dropbox (GaTech)/Asim/GMM Test/";
scg_path = '/home/dlin/Documents/MATLAB/GMM_SCG/';
addpath(genpath(scg_path))


files = dir(path + "DARPA nVNS Processed Output/*.mat");
for i = 24:length(files)
    file_i = files(i);
    sub = split(file_i.name, '_');
    sub = str2num(sub{1}(2:end));

    disp(string(datetime) + ": Processing " + sub)
    filename = "S" + num2str(sub) + "_processed_output_noLVET_SQIavail.mat";
    data = load(path + "DARPA nVNS Processed Output/" + filename);
    % Let's understand the file system
    findextrema = true;

    % I think hte first peak is hte rpeak time we'll also transpose to match
    % how I'm used to looking at signals
    scg_beats = data.scgBeats_chopped{1}(:, 2:end)';
    
    ao_range = [1, 250];
    ac_range = [250, 500];
    fs = 2000;
    disp("    finding AO")
    ao = sortPeaks(scg_beats, fs, 'axrange', ao_range, 'verbose', false, 'findextrema', findextrema);
    disp("    finding AC")
    [ac, tmp, score, nprob, d2] = sortPeaks(scg_beats, fs, 'axrange', ac_range, ...
        'nclusters', 3, 'verbose', false, 'findextrema', findextrema);
    ao(find(ao < 1)) = 1;
    % get most reasonable ao/ac 
    [~, best_ao_idx] = min(std(ao));
    [~, best_ac_idx] = min(std(ac));
    best_ao = ao(:, best_ao_idx);
    best_ac = ac(:, best_ac_idx);

    f = figure; 
    f.Position = [675 453 1500 750];
    subplot(3, 1, [1, 2])
    imagesc(normalize(scg_beats)); c= contrast(normalize(scg_beats)); colormap(c);
    hold on; p = plot(ao, 'r');  plot(ac, 'r');
    plot(best_ao, 'b'); plot(best_ac, 'b');

    ao_amps = general.indexMatrix(scg_beats, round(ao));
    ac_amps = general.indexMatrix(scg_beats, round(ac));
    subplot(3, 1, 3)
    beats = 1:100;
    p = plot(scg_beats(:, beats), 'k'); hold on;
    scatter(ao(beats, :), ao_amps(beats, :), 'r'); 
    scatter(ac(beats, :), ac_amps(beats, :), 'r')
    scatter(ao(beats, best_ao_idx), ao_amps(beats, best_ao_idx), 'b'); 
    scatter(ac(beats, best_ac_idx), ac_amps(beats, best_ac_idx), 'b')
    sgtitle("sortPeaks: Subject" + " " + sub)

    extrema_str = "";
    if findextrema == true
        extrema_str = "_findext";
    end
    saveas(f, path + "Figures/" + num2str(sub) + extrema_str + ".png")
    close all

    [maxes, mins] = general.slideTemplate(scg_beats, 30, 30, 2000, 'Axrange', ao_range, 'Closest', 'NumFeatures', 4);
    ao = [maxes;mins]';
    [maxes, mins] = general.slideTemplate(scg_beats, 30, 30, 2000, 'Axrange', ac_range, 'Closest', 'NumFeatures', 3);
    ac = [maxes;mins]';

        % get most reasonable ao/ac 
    [~, best_ao_idx] = min(std(ao));
    [~, best_ac_idx] = min(std(ac));
    best_ao = ao(:, best_ao_idx);
    best_ac = ac(:, best_ac_idx);

    f = figure; 
    f.Position = [675 453 1500 750];
    subplot(3, 1, [1, 2])
    imagesc(normalize(scg_beats)); c= contrast(normalize(scg_beats)); colormap(c);
    hold on; p = plot(ao, 'r');  plot(ac, 'r');
    plot(best_ao, 'b'); plot(best_ac, 'b');

    ao_amps = general.indexMatrix(scg_beats, round(ao));
    ac_amps = general.indexMatrix(scg_beats, round(ac));
    subplot(3, 1, 3)
    beats = 1:100;
    p = plot(scg_beats(:, beats), 'k'); hold on;
    scatter(ao(beats, :), ao_amps(beats, :), 'r'); 
    scatter(ac(beats, :), ac_amps(beats, :), 'r')
    scatter(ao(beats, best_ao_idx), ao_amps(beats, best_ao_idx), 'b'); 
    scatter(ac(beats, best_ac_idx), ac_amps(beats, best_ac_idx), 'b')
    sgtitle("slideTemplate: Subject" + " " + sub)

    saveas(f, path + "Figures/" + num2str(sub) + "_slide" + ".png")

    close all

end

% Notes
% 110 has crazy stuff happening beats 5425-5435
% 151 beat 4228 
%% want to iterate through the files


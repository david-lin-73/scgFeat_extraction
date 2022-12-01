function [maxlocs, minlocs, varargout] = slideTemplate(sigs, N_st, err, fs, varargin)
% -------------------------------------------------------------------------
% SUMMARY:
% This algorithm tracks multiple fiducial points using a sliding template.
% First the fiducial points are initialized using the initial tempalte
% (initial_template). Then for each beat a sliding template is made and its
% fiducial points are estimated using the previous fiducial points as
% references. From there the fiducial points of the beat being examined are
% estimated using the template associated with it. 
% ARGUMENTS (MANDATORY):
% - sigs             [NxM]     M signal segments of length N 
% - initial_template [Tx1]     template to estimate initial fid points
% - N_st             (dbl)     Number of beats to use in sliding template
% - err              (ms)      maximum range between old and new fid_point
% - fs                         sampling frequency
% ARGUMENTS (OPTIONAL)
% - Axrange         [1x2]      
% - To              (ms) Maximum fiducial point search range 
% - Closest         [flg]     Determine whether to select closest or max peak
% - NumFeatures               tracks ceil(NumFeatures/2) maxes anmd
%                                   1-ceil(NumFeatures/2) mins
% - EMD             [flg]     Used for better initialization of fid point
%                                   from template
% - Template        [flg]     Create a template using ens avg of first Nst
%                                   beats
% - Manual          [Bx2]     first column is index, second column is
%                                   whether it is a peak or valley (1/0)
% -
% OUTPUT:
% - maxlocs        [ceil(NumFeatures/2)x1]    Estimated fiducial point
%                                                   peaks from beats
% - minlocs        [1-ceil(NumFeatures/2)x1]  Estunated fudycuak point
%                                                   valleys from beats
%
% OUTPUT (OPTIONAL)
% - buff_templates [NxM]                      Sliding templates
% - temp_maxes     [ceil(NumFeatures/2)x1]    Estimated fiducial point
%                                                   peaks from templates
% - temp_mins      [1-ceil(NumFeatures/2)x1]  Estimated fiducial point
%                                                   valleys from templates
%
% USAGE:
% [maxlocs, minlocs, varargout] = slideTemplate(sigs, N_st, err, fs, varargin)
% -------------------------------------------------------------------------

if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp('NumFeatures', varargin{arg}); numfeats = varargin{arg+1}; end
        if strcmp('Emd', varargin{arg}); use_emd = 1; end
        if strcmp('Closest', varargin{arg}); closest = 1; end
        if strcmp('Axrange', varargin{arg}); axrange = varargin{arg+1}; axrange = round(axrange * fs/1000); end
        if strcmp('Template', varargin{arg}); initial_template = varargin{arg+1}; make_template = 0; end
        if strcmp('Nan', varargin{arg}); remove_nan = 1; end
        if strcmp('Remove', varargin{arg}); remove_range = round(varargin{arg+1}*fs/1000); end
        if strcmp('Manual', varargin{arg}); manual = 1; custom_idx = varargin{arg+1}; end
    end
end
if ~exist('numfeats', 'var'); numfeats = 2; end % number of features to track
if ~exist('use_emd', 'var'); use_emd = 0; end % emd for template enhancement 
if ~exist('closest', 'var'); closest = 1; end
if ~exist('make_template', 'var'); make_template = 1; end % make a template using ensemble averaging
if ~exist('remove_nan', 'var'); remove_nan = 0; end
if ~exist('axrange', 'var'); axrange = [1, size(sigs, 1)]; end
if ~exist('remove', 'var'); remove_range = 0; end
if ~exist('manual', 'var'); manual = 0; end


if make_template
    initial_template = sum(sigs(:, 1:N_st), 2);
end
% hold templates
buff_templates = zeros(size(sigs));


if ~manual
    
    % run emd to get a better idea of where the the fiducial points should be
    if use_emd
       [test, ~] = emd(initial_template, 'Display', 0);
       initial_template = test(:, 1);
    end

    % first detect the most prominant maxes and mins in template within the range specified
    [start_maxes, start_mins] = findMaxima(initial_template, numfeats, axrange, remove_range);
    
else

    start_maxes = custom_idx(find(custom_idx(:, 2) == 1), 1);
    start_mins = custom_idx(find(custom_idx(:, 2) == 0), 1);
    
end 


% place holders to track features
tot_maxes = zeros(length(start_maxes), size(sigs, 2));
tot_mins = zeros(length(start_mins), size(sigs, 2));

% place holders to track whether no peaks were found 
same_maxes = tot_maxes;
same_mins =  tot_mins;

% track the sliding templates fiducial points
temp_maxes = tot_maxes;
temp_mins = tot_mins;

% initialize the first fiducial point references to those found in template
old_sld_maxes = start_maxes;
old_sld_mins = start_mins;

% holds signals to create sliding template (initialize with 1st N_st sigs)
buff_signals = sigs(:, 1:N_st);

% initialization 
for beat = 1:N_st
    

    % ensemble average 
    new_sld_template = sum(buff_signals, 2)/N_st;


    buff_templates(:, beat) = new_sld_template;
    
    % first find fiducial point in new sliding template 
    new_sld_maxes = nearestPeak(old_sld_maxes, new_sld_template,...
        err, fs, 'Max', axrange, 'Closest', 1);
    new_sld_mins = nearestPeak(old_sld_mins, new_sld_template,...
        err, fs, 'Min', axrange, 'Closest', 1);

    % find fiducial point in examined beat
    tot_maxes(:, beat) = nearestPeak(new_sld_maxes, sigs(:, beat),...
        err, fs, 'Max', axrange, 'Closest', closest);
    tot_mins(:, beat) = nearestPeak(new_sld_mins, sigs(:, beat),...
        err, fs, 'Min', axrange, 'Closest', closest);

    % save template maxes and mins 
    temp_maxes(:, beat) = new_sld_maxes;
    temp_mins(:, beat) = new_sld_mins;

    % update the new AO/AC reference points
    old_sld_maxes = new_sld_maxes;
    old_sld_mins = new_sld_mins;
    
end

% initialize std and mean of fiducial points in template (used to identify
% when an incorrect peak is selected)
std_max = std(temp_maxes(:, 1:N_st),0,  2);
mean_max = mean(temp_maxes(:, 1:N_st), 2);
std_min = std(temp_mins(:, 1:N_st), 0, 2);
mean_min = mean(temp_mins(:, 1:N_st), 2);

for beat = N_st+1:size(sigs, 2)

        % create new template 
        new_sld_template = sum(buff_signals, 2)/N_st;
           

        % save template
        buff_templates(:, beat) = new_sld_template;


        % first find fiducial points in template
        new_sld_maxes = nearestPeak(old_sld_maxes, new_sld_template,...
            err, fs, 'Max', axrange, 'Closest', 1);%, 'Outlier', {mean_max, std_max});
        new_sld_mins = nearestPeak(old_sld_mins, new_sld_template,...
            err, fs, 'Min', axrange, 'Closest', 1);%, 'Outlier', {mean_min, std_min});

        % find fiducial points in examined beat and save whether no values
        % were found and the peak was just the reference
        [tot_maxes(:, beat), same_maxes(:, beat)] = nearestPeak(new_sld_maxes, sigs(:, beat),...
            err, fs, 'Max', axrange, 'Closest', closest);
        [tot_mins(:, beat), same_mins(:, beat)] = nearestPeak(new_sld_mins, sigs(:, beat),...
            err, fs, 'Min', axrange, 'Closest', closest);

        % save template maxes and mins 
        temp_maxes(:, beat) = new_sld_maxes;
        temp_mins(:, beat) = new_sld_mins;
        
        % update the new AO/AC reference points
        old_sld_maxes = new_sld_maxes;
        old_sld_mins = new_sld_mins;
        
        % update means and stds using the previous Nst beats
        std_max = std(temp_maxes(:, beat-N_st:beat),0,  2);
        mean_max = mean(temp_maxes(:, beat-N_st:beat), 2);
        std_min = std(temp_mins(:, beat-N_st:beat), 0, 2);
        mean_min = mean(temp_mins(:, beat-N_st:beat), 2);

        % update buffer by adding the signal just analyzed into buffer
        buff_signals = [buff_signals(:, 2:end), sigs(:, beat)];
        
end

if remove_nan
    for i = 1:size(tot_maxes, 1)
        if sum(same_maxes(i, :)) > 0
            tot_maxes;
        end
       tot_maxes(i, same_maxes(i, :) == 1) = nan;
    end
    
    for i = 1:size(tot_mins, 1)
       tot_mins(i, same_mins(i, :) == 1) = nan;
    end    
end

maxlocs = round(tot_maxes);        
minlocs = round(tot_mins);
varargout{1} = same_maxes; 
varargout{2} = same_mins;
varargout{3} = round(temp_maxes);
varargout{4} = round(temp_mins);
varargout{5} = buff_templates;
end

function [peakloc, varargout] = nearestPeak(fid_point, curr_seg, err, fs, findmax, axrange, varargin)
% -------------------------------------------------------------------------
% SUMMARY:
% This algorithm uses the old fiducial points (fid_point) as references to
% find the new fiducial point locations in the segment begin examined
% (curr_seg). If there are no candidate peaks/valleys the new fiducial
% point is set as the previous one. If there are candidate points, the
% fiducial point is chosen as the closest peak or the largest peak. Then
% chosen fiducial points are checked to make sure they aren't outliers.
%
% ARGUMENTS (MANDATORY):
% - fid_point       [Nx1]     fiducial points indices to look around
% - curr_seg        [Mx1]     segment being examined 
% - err             (ms)      maximum range between old and new fid_point
% - Fs              (dbl)     sampling frequency 
% - findmax         (str)     find maximum or minimum
% ARGUMENTS (OPTIONAL)
% - From            (samples) Minimum fiducial point search range 
% - To              (samples) Maximum fiducial point search range 
% - Closest         [flg]     Determine whether to select closest or max peak
% - Outlier
%                   {1}:    Mean vector of past fiducial points
%                   {2}:    Std vector of past fiducial points
% OUTPUT:
% - peak_loc        [Nx1]    Estimated fiducial points 
%
% USAGE:
% peakloc = nearestPeak(fid_point, curr_seg, err, 'Closest', 1)
% -------------------------------------------------------------------------
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp('Closest', varargin{arg}); closest = varargin{arg+1}; end
        if strcmp('Outlier', varargin{arg}); outlier = varargin{arg+1}; test_outlier = 1; end
    end
end
if ~exist('closest', 'var'); closest = 0; end
if ~exist('test_outlier', 'var'); test_outlier = 0; outlier = cell(2, 1); end

if strcmp('Max', findmax)
    findmax = 1;
else
    findmax = 0;
end

from = axrange(1);
to = axrange(2);

% extract the mean and std to detect outliers 
means = outlier{1};
stds = outlier{2};

% minimum std (for cases where the change in points is soo large that 
min_std = 1 * fs/1000;

% placeholder for peak locations 
peakloc = zeros(length(fid_point), 1);

% placeholder for whether reference value is chosen (initial guess will all
% not be forced to be the same)
sameloc = zeros(length(fid_point), 1);

% normalize segment
curr_seg = normalize(curr_seg);

% convert err from ms to samples
err = err * fs/1000;

% if looking for minimums just flip the signal to make it easier
if ~findmax
    curr_seg = -1*curr_seg;
end
[~, cand] = findpeaks(curr_seg);
cand(cand < from | cand > to) = [];

% for each point find the closest maxima in segment from fid_point 
for point = 1:length(fid_point)
        
    % first get the best candidates which is err away from the
    % fid_point
    best_cand = cand(abs(cand - fid_point(point)) < err); 
    
    % if no candidates are available just keep the fiducial point the same
    if isempty(best_cand)            

            peakloc(point) = fid_point(point);
            sameloc(point) = 1;

    % otherwise find the either the closest point or the point that gives
    % the largest amplitude
    else
         
        if closest
            
            % get the point(s) that are closest to the fiducial point
            [~, possible_points_idx] = min(abs(best_cand - fid_point(point)));
            possible_points = best_cand(possible_points_idx);
            
            % in case of ties find the maximum of the 2 to break ties
            peakloc(point) = possible_points(curr_seg(possible_points) ==...
                max(curr_seg(possible_points)));

        % or find largest max or the candidates within range
        else

            peakloc(point) = best_cand(curr_seg(best_cand) == max(curr_seg(best_cand)));
        
        end
        
        % if we are testing for outliers (if the new point is too far)
        if test_outlier 
            
            % if the selected point is 3x the std from the mean of previous
            % points correct the point 
            if abs(peakloc(point) - means(point)) > max(3*stds(point), 3*min_std)
                scalar = sign(peakloc(point)-means(point))*0.1;
                peakloc(point) = fid_point(point) + scalar*max(std(point), min_std);
            end
            
        end

    end
 
end

varargout{1} = sameloc;
end

% This is def something 
function [maxes, mins] = findMaxima(signal, numfeats, axrange, remove_range)
% find the extrema given a signal 

% first normalize signal 
signal = normalize(signal);

% calculate the number of features to find within range
nummaxes = ceil(numfeats/2); nummins = numfeats - nummaxes;
maxes = zeros(nummaxes, 1); mins = zeros(nummins, 1);


% find maximas and minimas 
[mvals, mlocs] = findpeaks(signal, 'MinPeakDistance', remove_range);
[~, vlocs] = findpeaks(-1*signal, 'MinPeakDistance', remove_range);
vvals = signal(vlocs);

% remove locations 
poor_mlocs = mlocs > axrange(2) | mlocs < axrange(1);
poor_vlocs = vlocs > axrange(2) | vlocs < axrange(1);
mlocs(poor_mlocs) = []; 
vlocs(poor_vlocs) = [];
mvals(poor_mlocs) = []; 
vvals(poor_vlocs) = [];
    
% use this if you want to choose points that are considered maximas 
[~, mrank] = sort(mvals, 'descend'); [~, vrank] = sort(vvals, 'ascend');
if ~isempty(mvals); maxes = mlocs(mrank(1:min(nummaxes, length(mvals)))); end
if ~isempty(vvals); mins = vlocs(vrank(1:min(nummins, length(vvals)))); end


end
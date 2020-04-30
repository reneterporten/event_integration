function [resulting_data] = rt_searchlight(cfgdata, data)

% data should be in format:
% data = {erp_A_pre; erp_B_pre; erp_C_pre};

timewidth       = cfgdata.timewidth;
timesteps       = cfgdata.timesteps;
neighbours      = cfgdata.neighbours;
searchspace     = cfgdata.searchspace;
searchtime      = cfgdata.searchtime;
avgovertrials   = cfgdata.avgovertrials;
avgovertime     = cfgdata.avgovertime;
avgoverspace    = cfgdata.avgoverspace;

% Default check values for the data. Will be true if requested timestep do
% not match temporal frequency of the data
timestep_warning = false;

% Calculate timewindow that can be fitted to the actual data
min_time    = min(data{1,1}.time);
max_time    = max(data{1,1}.time);
steps_time  = data{1,1}.time(2)-data{1,1}.time(1);

if timesteps < steps_time
    warning('Indicated time steps of shifting window cannot be fitted to actual time steps of data. Timesteps will be approximated closest to indicated value.')
    timestep_warning = true;
elseif ~rem(timesteps,steps_time)*timesteps/steps_time == 0
    warning('Indicated time steps of shifting window cannot be fitted to actual time steps of data. Timesteps will be approximated closest to indicated value.')
    timestep_warning = true;
end

% Approximate time steps of window 
req_timesteps       = min_time + timesteps;
[val,idx]           = min(abs(data{1,1}.time-req_timesteps));
use_timestepsidx    = idx - 1; % Index by which the timewindow needs to shift to match requested timestep  

% Approximate time width of window 
req_timewidth       = min_time + timewidth;
[val,idx]           = min(abs(data{1,1}.time-req_timewidth));
use_timewidthidx    = idx - 1; % Index by which the timewindow needs to shift to match requested timestep  

% Calculate number of time windows that fit the data
% These steps can be adjusted such that the fit of the last timewindow is
% set to some threshold. (e.g. include timewindow if 90% of data is preserved)
window_counter = 0;
for window_fit_vec = 1:use_timestepsidx:length(data{1,1}.time) 
    if (window_fit_vec + use_timewidthidx) <= length(data{1,1}.time)
        window_counter = window_counter + 1;
    end  
end

% Store number of trials for each condition
data_trial = zeros(size(data, 1),1);
for dataNr = 1:size(data_trial,1)
    data_trial(dataNr,1) = size(data{dataNr,1}.trialinfo, 1);
end

% Loop through channels and extract relevant data over time
% Preallocate final structure
resulting_data                  = [];
resulting_data.searchlight      = zeros(length(data{1,1}.label), window_counter, sum(data_trial, 1), sum(data_trial, 1));
for idxLabel = 1:length(data{1,1}.label)
    
    disp('Extracting channel indices.')
    indicesChannel = [idxLabel];
    % For each channel label check for neighbours and store their indices
    for neighb = 1:length(neighbours(idxLabel).neighblabel)
        indicesChannel = [indicesChannel, find(strcmp(data{1,1}.label, neighbours(idxLabel).neighblabel(neighb)))];
    end
    
    time_start_idx = 1;
    data_search = [];
    time_search = [];
    disp('Extracting data per timewindow.')
    % For each channel and its including neighbours loop through time
    % windows
    for timewindow_nr = 1:1%window_counter
        
        disp(strcat('Calculating timewindow:', int2str(timewindow_nr)))
        
        % Create Trial structure including all conditions
        startTrial = 0;
        for condNr = 1:size(data_trial, 1)
            for nrTrial = 1:data_trial(condNr)
                data_search(startTrial + nrTrial, [1:length(indicesChannel)], 1:(use_timewidthidx + 1)) = data{condNr}.trial(nrTrial, indicesChannel,time_start_idx:(time_start_idx+use_timewidthidx));
                time_search(timewindow_nr, 1:(use_timewidthidx + 1)) = data{condNr}.time(time_start_idx:(time_start_idx+use_timewidthidx)); 
            end
            startTrial = startTrial + data_trial(condNr);
        end

        % Do comparison of trial x trials
        % The trial structure is derived from the rows of data_search
        % data_search includes trial x channels(incl neighbours) x time
        temp_comp = zeros(size(data_search, 1), size(data_search, 1));
        for rowData = 1:size(data_search, 1)
            for colData = 1:size(data_search, 1)
                % This results in a symmetrical trial x trial matrix for
                % this particular channel & timewindow combination
                temp_comp(rowData, colData) = corr2(squeeze(data_search(rowData, :, :)), squeeze(data_search(colData, :, :)));
            end   
        end
        
        % Store temporary data of trial comparison to matrix
        resulting_data.searchlight(idxLabel, timewindow_nr, :, :) = temp_comp;
        
        time_start_idx = time_start_idx + use_timestepsidx;
    end
        
    disp(strcat('Calculating label:', int2str(idxLabel)))
    
end

%disp('Store data into structure.')
resulting_data.windownr_time    = time_search;
resulting_data.label            = data{1,1}.label;
resulting_data.dimord           = 'chan_windownr_trial_trial';





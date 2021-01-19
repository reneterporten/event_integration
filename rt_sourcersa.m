function [dataRSA] = rt_sourcersa(cfgdata, data)

time            = cfgdata.time;
atlas           = cfgdata.atlas;
parcelrm        = cfgdata.parcelrm;
timewidth       = cfgdata.timewidth;
timesteps       = cfgdata.timesteps;
searchspace     = cfgdata.searchspace;
searchtime      = cfgdata.searchtime;
avgovertrials   = cfgdata.avgovertrials;
avgovertime     = cfgdata.avgovertime;
avgoverspace    = cfgdata.avgoverspace;

% Default check values for the data. Will be true if requested timestep do
% not match temporal frequency of the data
timestep_warning = false;

% Calculate timewindow that can be fitted to the actual data
min_time    = min(time);
max_time    = max(time);
steps_time  = time(2)-time(1);

if timesteps < steps_time
    warning('Indicated time steps of shifting window cannot be fitted to actual time steps of data. Timesteps will be approximated closest to indicated value.')
    timestep_warning = true;
elseif ~rem(timesteps,steps_time)*timesteps/steps_time == 0
    warning('Indicated time steps of shifting window cannot be fitted to actual time steps of data. Timesteps will be approximated closest to indicated value.')
    timestep_warning = true;
end

% Approximate time steps of window 
req_timesteps       = min_time + timesteps;
[val,idx]           = min(abs(time-req_timesteps));
use_timestepsidx    = idx - 1; % Index by which the timewindow needs to shift to match requested timestep  

% Approximate time width of window 
req_timewidth       = min_time + timewidth;
[val,idx]           = min(abs(time-req_timewidth));
use_timewidthidx    = idx - 1; % Index by which the timewindow needs to shift to match requested timestep  

% Calculate number of time windows that fit the data
% These steps can be adjusted such that the fit of the last timewindow is
% set to some threshold. (e.g. include timewindow if 90% of data is preserved)
window_counter = 0;
for window_fit_vec = 1:use_timestepsidx:length(time) 
    if (window_fit_vec + use_timewidthidx) <= length(time)
        window_counter = window_counter + 1;
    end  
end

% Get all unique parcels and skip the 0 entries
allParcels = unique(atlas.parcellation1D);
allParcels = allParcels(2:end);

% Remove parcels
if ~isempty(parcelrm)
    allParcels = allParcels(~ismember(allParcels, parcelrm));
end

resulting_data                  = [];
resulting_data.searchlight      = zeros(length(allParcels), window_counter, length(data), length(data));
% Initiate loop to go through all parcels and time-windows
for idxParc = 1:length(allParcels)
    
    % Extract indices of parcel
    indicesParcel = atlas.parcellation1D == allParcels(idxParc);
    
    % Loop through each time-window
    time_start_idx  = 1;
    data_search     = [];
    for timewindow_nr = 1:window_counter
        
        % For each condition select the parcels and timewindows
        for condNr = 1:size(data, 1)
            % Parcel specific data
            parcel_data = data{condNr}(indicesParcel,:);
            for parcels = 1:size(parcel_data,1)
                data_search(condNr, parcels, 1:(use_timewidthidx + 1)) = parcel_data(parcels, time_start_idx:(time_start_idx+use_timewidthidx));
            end
        end
        
        % Advance in time window
        time_start_idx = time_start_idx + use_timestepsidx;
        
        % Correlate neural data against conditions, results in neural RSM
        data_search_2d  = data_search(:,:);
        temp_comp       = corr(data_search_2d');
        
        % For now store the data
        resulting_data.searchlight(idxParc, timewindow_nr, :, :) = temp_comp;
        
    end
    
    disp(strcat('Calculating label:', int2str(idxParc)))
    
end

% Compare neural RSM to prediction rsm
sanitydata      = xlsread('/project/3012026.13/scripts_RT/sanity_pre.xlsx'); % Visual similarity
% Create prediction RDM that matches trial x trial structure of data
cfg             = [];
cfg.offdiag     = true; % Use postdata structure for off diagional (phase comparison)
predictionRDM   = rt_predictionRDMxlsx(cfg, sanitydata, sanitydata);
predictionRDM(isnan(predictionRDM)) = -1; % Use -1 instead of NaN for visual similarity

predictionRDM_vec   = vectorizeSimmat(predictionRDM);
dataRSA             = [];
dataRSA.RSM         = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.parcels     = allParcels;

tic
for nrChan = 1:size(resulting_data.searchlight, 1)
    disp(strcat('Calculating Channel:', int2str(nrChan)))
    for nrTimewin = 1:size(resulting_data.searchlight, 2)

        currentNeural               = squeeze(resulting_data.searchlight(nrChan, nrTimewin, :, :));
        currentNeural               = vectorizeSimmat(currentNeural);
        % dataRSA is the resulting structure containing chan x timewindow
        % information
        dataRSA.RSM(nrChan, nrTimewin)  = corr(predictionRDM_vec', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');

    end
end
toc







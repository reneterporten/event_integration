function rt_searchlight(cfgdata, data, subj)

disp(strcat('Searchlight on subject: ', subj))

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
resulting_data.searchlight      = zeros(length(neighbours), window_counter, sum(data_trial, 1), sum(data_trial, 1));
resulting_data_GLM              = [];
resulting_data_GLM.searchlight  = zeros(length(neighbours), window_counter, 6, 6);
tic
for idxLabel = 1:length(neighbours) % This loops through 271 labels
    
    %disp('Extracting channel indices.')
    indicesChannel_dH = [idxLabel];
    indicesChannel_dV = [find(strcmp(data{1,1}.label, {strcat(neighbours(idxLabel).label, '_dV')}))];
    % For each channel label check for neighbours and store their indices
    for neighb = 1:length(neighbours(idxLabel).neighblabel)
        % Indicate channel
        neigh_chan_dH       = {strcat(neighbours(idxLabel).neighblabel{neighb}, '_dH')};
        neigh_chan_dV       = {strcat(neighbours(idxLabel).neighblabel{neighb}, '_dV')};
        indicesChannel_dH   = [indicesChannel_dH, find(strcmp(data{1,1}.label, neigh_chan_dH))];
        indicesChannel_dV   = [indicesChannel_dV, find(strcmp(data{1,1}.label, neigh_chan_dV))];
    end
    
    % Combine the horizontal and vertical channel indices
    indicesChannel = [indicesChannel_dH, indicesChannel_dV];
    
    time_start_idx = 1;
    data_search = [];
    time_search = [];
    %disp('Extracting data per timewindow.')
    % For each channel and its including neighbours loop through time
    % windows
    for timewindow_nr = 1:window_counter
        
        %disp(strcat('Calculating timewindow:', int2str(timewindow_nr)))
        
        % Create Trial structure including all conditions
        %tic
        startTrial = 0;
        for condNr = 1:size(data_trial, 1)
            for nrTrial = 1:data_trial(condNr)
                data_search(startTrial + nrTrial, [1:length(indicesChannel)], 1:(use_timewidthidx + 1)) = data{condNr}.trial(nrTrial, indicesChannel,time_start_idx:(time_start_idx+use_timewidthidx));
                time_search(timewindow_nr, 1:(use_timewidthidx + 1)) = data{condNr}.time(time_start_idx:(time_start_idx+use_timewidthidx)); 
            end
            startTrial = startTrial + data_trial(condNr);
        end
        %stopWatch = toc;
        %disp(strcat('Time to create trial structure: ', num2str(stopWatch)))
        
        % Subtract the channel specific mean from the data
        data_search = rt_subtrmean(data_search);

        % Do comparison of trial x trials
        % The trial structure is derived from the rows of data_search
        % data_search includes trial x channels(incl neighbours) x time
        data_search_2d  = data_search(:,:);
        temp_comp       = corr(data_search_2d');
        
        % Store temporary data of trial comparison to matrix
        resulting_data.searchlight(idxLabel, timewindow_nr, :, :) = temp_comp;
        
        % Store data for GLM calculation (data extracted)
        resulting_data_GLM.searchlight(idxLabel, timewindow_nr, :, :) = rt_extractData(temp_comp);
        
        time_start_idx = time_start_idx + use_timestepsidx;
        
    end
      
    disp(strcat('Calculating label:', int2str(idxLabel)))
    
end
toc

% Save extracted data for GLM
mkdir(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj))
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'resulting_data_GLM.mat'), 'resulting_data_GLM', '-v7.3');

% Compare neural data to prediction RSMs

% Visual similarity 1
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'visual_similarity_1.mat'), 'visual_similarity_1')
predictionRDM_vec_vis1       = vectorizeSimmat(visual_similarity_1);
clear visual_similarity_1
% Visual similarity 2
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_vs2_magn.mat']), 'vs2_magn')
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_vs2_phas.mat']), 'vs2_phas')
predictionRDM_vec_magn       = vectorizeSimmat(vs2_magn);
predictionRDM_vec_phas       = vectorizeSimmat(vs2_phas);
clear vs2_magn
clear vs2_phas
% Narrative insight
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'narrativeInsight.mat'), 'narrativeInsight')
predictionRDM_vec_narinsight = vectorizeSimmat(narrativeInsight);

% Create narrative insight pred RSM with medium & high noise levels
% noise_noeffect    = randi([600 900], 432, 432);
% noise_noeffect    = mednoise./1000;
% effectdir           = randi(2,432)-1;
% effectdir(~effectdir)=-1;
% noise_effect   = randi([300 600], 432, 432);
% noise_effect   = highnoise./1000;
% predictionRDM_vec_narinsight_noise_effect   = vectorizeSimmat(narrativeInsight.*noise_effect);
% predictionRDM_vec_narinsight_noise_noeffect  = vectorizeSimmat((narrativeInsight.*effectdir).*noise_noeffect);

% Create narrative insight pred RSM with only the pot phase expectations
predictionRDM_vec_narinsight_prepost = narrativeInsight(217:432,217:432);
predictionRDM_vec_narinsight_prepost = vectorizeSimmat(predictionRDM_vec_narinsight_prepost);

clear narrativeInsight
% Behavior insight
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_behaviorInsight.mat']), 'behaviorInsight')
predictionRDM_vec_behavinsight = vectorizeSimmat(behaviorInsight);
clear behaviorInsight
% Narrative mismatch
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_narrativeMismatch.mat']), 'narrativeMismatch')
predictionRDM_vec_narmism = vectorizeSimmat(narrativeMismatch);
clear narrativeMismatch
% (Un)linking mechanism
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'un_linkMech.mat'), 'un_linkMech')
predictionRDM_vec_unlink = vectorizeSimmat(un_linkMech);
clear un_linkMech

% Prepare dataRSA structures
dataRSA                 = [];
% Subfields for original structure
dataRSA.vis1            = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.magn            = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.phas            = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.narinsight      = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.behavinsight    = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.narmism         = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.unlink          = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));

% Subfields for noise structures
% dataRSA.narinsight_mednoise    = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
% dataRSA.narinsight_highnoise   = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));

% Subfields for pre and post phase splitting;
dataRSA.narinsight_pre     = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.narinsight_post    = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));

% Subfields for permuted structure
% dataRSA.vis1_rand           = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
% dataRSA.magn_rand           = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
% dataRSA.phas_rand           = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
% dataRSA.narinsight_rand     = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
% dataRSA.behavinsight_rand   = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
% dataRSA.narmism_rand        = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
% dataRSA.unlink_rand         = zeros(size(resulting_data.searchlight,1), size(resulting_data.searchlight,2));
dataRSA.windownr_time       = time_search;
dataRSA.dimord              = 'chan_windownr';
dataRSA.label               = data{1,1}.label;

% Create permuted structure
% disp('Create permuted structure')
% num_perm = 500;
% n_trials = squeeze(resulting_data.searchlight(1, 1, :, :));
% n_trials = length(vectorizeSimmat(n_trials));
% perm_idx = zeros(num_perm, n_trials);
% rng(5)
% for r = 1:num_perm
%     perm_idx(r,:) = randperm(n_trials);
% end

tic
for nrChan = 1:size(resulting_data.searchlight, 1)
    disp(strcat('Calculating Channel:', int2str(nrChan)))
    for nrTimewin = 1:size(resulting_data.searchlight, 2)
        currentNeural               = squeeze(resulting_data.searchlight(nrChan, nrTimewin, :, :));
        currentNeural_pre           = currentNeural(1:216,1:216); % extract only pre data
        currentNeural_post          = currentNeural(217:432,217:432); % extract only pre data
        currentNeural               = vectorizeSimmat(currentNeural);
        currentNeural_pre           = vectorizeSimmat(currentNeural_pre);
        currentNeural_post          = vectorizeSimmat(currentNeural_post);
        
        % dataRSA is the resulting structure containing chan x timewindow
        % information
        disp('Correlating model x actual data')
        dataRSA.vis1(nrChan, nrTimewin)         = corr(predictionRDM_vec_vis1', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
        dataRSA.magn(nrChan, nrTimewin)         = corr(predictionRDM_vec_magn', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
        dataRSA.phas(nrChan, nrTimewin)         = corr(predictionRDM_vec_phas', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
        dataRSA.narinsight(nrChan, nrTimewin)   = corr(predictionRDM_vec_narinsight', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
        dataRSA.behavinsight(nrChan, nrTimewin) = corr(predictionRDM_vec_behavinsight', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
        dataRSA.narmism(nrChan, nrTimewin)      = corr(predictionRDM_vec_narmism', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
        dataRSA.unlink(nrChan, nrTimewin)       = corr(predictionRDM_vec_unlink', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
        
        % Correlate with noisy interaction matrix
%         dataRSA.narinsight_mednoise(nrChan, nrTimewin)  = corr(predictionRDM_vec_narinsight_mednoise', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
%         dataRSA.narinsight_highnoise(nrChan, nrTimewin) = corr(predictionRDM_vec_narinsight_highnoise', currentNeural', 'Type', 'Spearman', 'Rows', 'complete');
%         
        % Correlate post phase expectations with splitted data
        dataRSA.narinsight_pre(nrChan, nrTimewin)  = corr(predictionRDM_vec_narinsight_prepost', currentNeural_pre', 'Type', 'Spearman', 'Rows', 'complete');
        dataRSA.narinsight_post(nrChan, nrTimewin) = corr(predictionRDM_vec_narinsight_prepost', currentNeural_post', 'Type', 'Spearman', 'Rows', 'complete');
        
%         % Calculate correlation between model and permuted data
%         disp('Permute neural data and average')
%         currentNeural_rand = zeros(size(currentNeural));
%         for randNr = 1:num_perm
%             currentNeural_rand = currentNeural_rand + currentNeural(perm_idx(randNr,:));
%         end
%         currentNeural_rand = currentNeural_rand/num_perm;
%         
%         % Data RSA structure from the model and permuted comparisons
%         disp('Correlating model x permuted data')
%         dataRSA.vis1_rand(nrChan, nrTimewin)         = corr(predictionRDM_vec_vis1', currentNeural_rand', 'Type', 'Spearman', 'Rows', 'complete');
%         dataRSA.magn_rand(nrChan, nrTimewin)         = corr(predictionRDM_vec_magn', currentNeural_rand', 'Type', 'Spearman', 'Rows', 'complete');
%         dataRSA.phas_rand(nrChan, nrTimewin)         = corr(predictionRDM_vec_phas', currentNeural_rand', 'Type', 'Spearman', 'Rows', 'complete');
%         dataRSA.narinsight_rand(nrChan, nrTimewin)   = corr(predictionRDM_vec_narinsight', currentNeural_rand', 'Type', 'Spearman', 'Rows', 'complete');
%         dataRSA.behavinsight_rand(nrChan, nrTimewin) = corr(predictionRDM_vec_behavinsight', currentNeural_rand', 'Type', 'Spearman', 'Rows', 'complete');
%         dataRSA.narmism_rand(nrChan, nrTimewin)      = corr(predictionRDM_vec_narmism', currentNeural_rand', 'Type', 'Spearman', 'Rows', 'complete');
%         dataRSA.unlink_rand(nrChan, nrTimewin)       = corr(predictionRDM_vec_unlink', currentNeural_rand', 'Type', 'Spearman', 'Rows', 'complete');
%         
    end
end
toc

mkdir(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj))
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'dataRSA.mat'), 'dataRSA', '-v7.3');














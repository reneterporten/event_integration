%% Apply RSA with spatial and temporal searchlight to timelock data
%% Default path

addpath /home/common/matlab/fieldtrip/
addpath /home/common/matlab/fieldtrip/qsub/
addpath /project/3012026.13/scripts_RT/
ft_defaults

root_dir     = '/project/3012026.13/processed/';
save_dir     = '/project/3012026.13/processed_RT/time_lock/';

subjects = {'sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-022','sub-023',...
    'sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-031'...
    'sub-032','sub-033','sub-034','sub-035','sub-036','sub-037','sub-038'};


%% Load in timelock data and apply spatial/temporal searchlight

for iSubject = 1:1%length(subjects)
    
    % Load Data
    subj                = subjects{iSubject};
    load(fullfile(save_dir, [subj,'_ERP_basecon']), 'erp_A_pre', 'erp_B_pre', 'erp_X_pre');
    my_data             = {erp_A_pre; erp_B_pre; erp_X_pre};
    
    % Get neighbour structure, do it only once
    if iSubject == 1
        cfg             = [];
        cfg.method      = 'triangulation';
        cfg.layout      = 'CTF275_helmet.mat';
        cfg.channel     = my_data{1,1}.label;
        neighbours      = ft_prepare_neighbours(cfg);
    end
    
    % Apply searchlight function
    cfg                 = [];
    cfg.timewidth       = 0.200; % Width of shifting timewindow
    cfg.timesteps       = 0.050; % Steps of shifting timewindow
    cfg.neighbours      = neighbours;
    cfg.searchspace     = 'yes';
    cfg.searchtime      = 'yes';
    cfg.avgovertrials   = 'yes';
    cfg.avgovertime     = 'no';
    cfg.avgoverspace    = 'no';
    sl_data_pre         = rt_searchlight(cfg, my_data);
       
end

keep iSubject subjects root_dir save_dir cfgdata sl_data_pre my_data


%% Create prediction RDM and compare to neural RDM

my_template     = [NaN, -1, 1; -1, NaN, 1; 1, 1, NaN];
predictionRDM   = rt_predictionRDM(my_template, my_data);
neuralRDM       = sl_data_pre.searchlight;

predictionRDM_vec   = vectorizeSimmat(predictionRDM);
dataRSA             = zeros(size(neuralRDM,1), size(neuralRDM,2));
for nrChan = 1:size(neuralRDM, 1)
    disp(strcat('Calculating Channel:', int2str(nrChan)))
    for nrTimewin = 1:size(neuralRDM, 2)

        currentNeural               = squeeze(neuralRDM(nrChan, nrTimewin, :, :));
        currentNeural               = vectorizeSimmat(currentNeural); 
        dataRSA(nrChan, nrTimewin)  = corr(predictionRDM_vec', currentNeural', 'Type', 'Kendall', 'Rows', 'complete');

    end
end









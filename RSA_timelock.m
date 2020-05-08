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

for iSubject = 1:length(subjects)
    
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
    cfg.timewidth       = 0.750; % Width of shifting timewindow
    cfg.timesteps       = 0.400; % Steps of shifting timewindow
    cfg.neighbours      = neighbours;
    cfg.searchspace     = 'yes'; % Not included in rt_searchlight yet
    cfg.searchtime      = 'yes'; % Not included in rt_searchlight yet
    cfg.avgovertrials   = 'yes'; % Not included in rt_searchlight yet
    cfg.avgovertime     = 'no'; % Not included in rt_searchlight yet
    cfg.avgoverspace    = 'no'; % Not included in rt_searchlight yet
    sl_data_pre         = rt_searchlight(cfg, my_data);
       
end

keep iSubject subjects root_dir save_dir cfgdata sl_data_pre my_data data

save(fullfile('/project/3012026.13/processed_RT/', 'searchlight_test2.mat'), 'sl_data', '-v7.3')


%% Create prediction RDM and compare to neural RDM

predata     = xlsread('/project/3012026.13/scripts_RT/prediction_pre_noBX.xlsx');
postdata    = xlsread('/project/3012026.13/scripts_RT/prediction_post_noBX.xlsx');
sanitydata  = xlsread('/project/3012026.13/scripts_RT/sanity_pre.xlsx');
% Create prediction RDM that matches trial x trial structure of data
predictionRDM   = rt_predictionRDMxlsx(predata, postdata);
neuralRDM       = sl_data.searchlight;

predictionRDM_vec   = vectorizeSimmat(predictionRDM);
dataRSA             = zeros(size(neuralRDM,1), size(neuralRDM,2));
for nrChan = 1:size(neuralRDM, 1)
    disp(strcat('Calculating Channel:', int2str(nrChan)))
    for nrTimewin = 1:size(neuralRDM, 2)

        currentNeural               = squeeze(neuralRDM(nrChan, nrTimewin, :, :));
        currentNeural               = vectorizeSimmat(currentNeural);
        % dataRSA is the resulting structure containing chan x timewindow
        % information
        dataRSA(nrChan, nrTimewin)  = corr(predictionRDM_vec', currentNeural', 'Type', 'Kendall', 'Rows', 'complete');

    end
end

% Integrate data from RSA into structure that is intepretable by fieldtrip
cfg                 = [];
cfg.avgoverrpt      = 'yes';
cfg.avgovertime     = 'yes';
sl_data_test        = ft_selectdata(cfg, my_data{1});
sl_data_test.trial  = dataRSA(:,4);

% Plot the data onto topography
cfg                 = [];
cfg.zlim            = [-0.06 0.06]; 
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'trial';
cfg.comment         = 'no';
figure; ft_topoplotER(cfg,sl_data_test); colorbar
set(gcf,'color','w')

ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu'))) 

% Plot channel x timewindow data
imagesc(dataRSA, [-0.06 0.06])
set(gcf,'color','w')
colormap(flipud(brewermap(64,'RdBu'))) 





load(fullfile('/project/3012026.13/processed_RT/', 'searchlight_test.mat'))


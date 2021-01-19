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


%% Load in all RSA data per subject

% Integrate data from RSA into structure that is intepretable by fieldtrip
load(fullfile('/project/3012026.13/processed_RT/', 'template_timelock.mat'));
cfg                 = [];
cfg.avgoverrpt      = 'yes';
cfg.avgovertime     = 'no';
template_data       = ft_selectdata(cfg, template_timelock);
template_data.time  = [1:15]*0.125;
template_data       = rmfield(template_data, 'trial');
template_data.cfg   = [];
clear template_timelock

% Create data structure that contains RSA data of all subjects
% Allocate overall structure
ignSub = [4, 10, 12, 20, 30, 32];
allDataRSA = cell((length(subjects)-length(ignSub)),1);

rsa_count = 1;
for subData = 1:length(subjects)
    
    if ~ismember(subData, ignSub)

        disp(strcat('Creating RSA structure:', int2str(subData)))

        subj = subjects{subData};

        load(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'dataRSA.mat'));

        data_RSA                 = template_data;
        data_RSA.vis1            = dataRSA.vis1;
        data_RSA.magn            = dataRSA.magn;
        data_RSA.phas            = dataRSA.phas;
        data_RSA.narinsight      = dataRSA.narinsight;
        data_RSA.behavinsight    = dataRSA.behavinsight;
        data_RSA.narmism         = dataRSA.narmism;
        data_RSA.unlink          = dataRSA.unlink;
        data_RSA.narinsight_pre   = dataRSA.narinsight_pre;
        data_RSA.narinsight_post  = dataRSA.narinsight_post;
        data_RSA.narinsight_mednoise  = dataRSA.narinsight_mednoise;
        data_RSA.narinsight_highnoise  = dataRSA.narinsight_highnoise;
        %data_RSA.windownr_time   = dataRSA.windownr_time;
        data_RSA.label           = template_data.label;

        allDataRSA{rsa_count}   = data_RSA;
        rsa_count               = rsa_count + 1;
    
    end

end

keep allDataRSA ignSub subjects save_dir root_dir


%% Calculate average scores for each prediction RSM

cfg                 = [];
disp('Vis1 avg...')
cfg.parameter       = 'vis1';
vis1_avg            = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('Magn avg...')
cfg.parameter       = 'magn';
magn_avg            = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('Phas avg...')
cfg.parameter       = 'phas';
phas_avg            = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('Narinsight avg...')
cfg.parameter       = 'narinsight';
narinsight_avg      = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('Behavinsight avg...')
cfg.parameter       = 'behavinsight';
behavinsight_avg    = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('Narmism avg...')
cfg.parameter       = 'narmism';
narmism_avg         = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('Unlink avg...')
cfg.parameter       = 'unlink';
unlink_avg          = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('narinsight_pre avg...')
cfg.parameter       = 'narinsight_pre';
narinsight_pre_avg          = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('narinsight_post avg...')
cfg.parameter       = 'narinsight_post';
narinsight_post_avg          = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('narinsight_mednoise avg...')
cfg.parameter       = 'narinsight_mednoise';
narinsight_mednoise_avg          = ft_timelockgrandaverage(cfg, allDataRSA{:});
disp('narinsight_highnoise avg...')
cfg.parameter       = 'narinsight_highnoise';
narinsight_highnoise_avg          = ft_timelockgrandaverage(cfg, allDataRSA{:});

disp('Saving avg structures...')
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'vis1_avg.mat'), 'vis1_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'magn_avg.mat'), 'magn_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'phas_avg.mat'), 'phas_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'narinsight_avg.mat'), 'narinsight_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'behavinsight_avg.mat'), 'behavinsight_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'narmism_avg.mat'), 'narmism_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'unlink_avg.mat'), 'unlink_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'narinsight_pre_avg.mat'), 'narinsight_pre_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'narinsight_post_avg.mat'), 'narinsight_post_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'narinsight_mednoise_avg.mat'), 'narinsight_mednoise_avg', '-v7.3');
save(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'narinsight_highnoise_avg.mat'), 'narinsight_highnoise_avg', '-v7.3');


%% Load average data

load(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'vis1_avg.mat'));
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'magn_avg.mat'));
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'phas_avg.mat'));
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'narinsight_avg.mat'));
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'behavinsight_avg.mat'));
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'narmism_avg.mat'));
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', 'unlink_avg.mat'));


%% Plot results

% Topographies
% Plot the data onto topography
cfg                 = [];
cfg.xlim            = [1:15]*0.125;
cfg.zlim            = [-0.01 0.01];
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'avg';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, vis1_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, magn_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, phas_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, narinsight_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, behavinsight_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, narmism_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, unlink_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, narinsight_pre_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, narinsight_post_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))



figure; ft_topoplotER(cfg, narinsight_mednoise_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, narinsight_highnoise_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

% ------


% Plot channel x timewindow data
figure('Renderer', 'painters', 'Position', [10 10 300 500])
imagesc(vis1_avg.avg, [-0.01 0.01])
xlabel('Timewindow')
ylabel('Channel')
title('Similarity (r)')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 

figure('Renderer', 'painters', 'Position', [10 10 300 500])
imagesc(magn_avg.avg)
xlabel('Timewindow')
ylabel('Channel')
title('Similarity (r)')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 

figure('Renderer', 'painters', 'Position', [10 10 300 500])
imagesc(phas_avg.avg)
xlabel('Timewindow')
ylabel('Channel')
title('Similarity (r)')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 

figure('Renderer', 'painters', 'Position', [10 10 300 500])
imagesc(narinsight_avg.avg)
xlabel('Timewindow')
ylabel('Channel')
title('Similarity (r)')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 

figure('Renderer', 'painters', 'Position', [10 10 300 500])
imagesc(behavinsight_avg.avg)
xlabel('Timewindow')
ylabel('Channel')
title('Similarity (r)')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 

figure('Renderer', 'painters', 'Position', [10 10 300 500])
imagesc(narmism_avg.avg)
xlabel('Timewindow')
ylabel('Channel')
title('Similarity (r)')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 

figure('Renderer', 'painters', 'Position', [10 10 300 500])
imagesc(unlink_avg.avg)
xlabel('Timewindow')
ylabel('Channel')
title('Similarity (r)')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 


%% Plot topos with mask

rangetop                = (length(narinsight_avg.avg(:))/100)*10;

narinsight_avg_top5     = narinsight_avg.avg;
[~,idx]                 = sort(narinsight_avg_top5(:), 'descend');
narinsight_toppos       = narinsight_avg_top5(idx(1:round(rangetop)));
[~,idx]                 = sort(narinsight_avg_top5(:), 'ascend');
narinsight_topneg       = narinsight_avg_top5(idx(1:round(rangetop)));

narinsight_masked = narinsight_avg;
topmask = ismember(narinsight_masked.avg, narinsight_toppos) + ismember(narinsight_masked.avg, narinsight_topneg);
topmask = topmask > 0;
narinsight_masked.avg(topmask) = 0;




% Topographies
% Plot the data onto topography
cfg                 = [];
cfg.xlim            = [1:15]*0.125;
cfg.zlim            = [-0.005 0.005];
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'avg';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, narinsight_masked); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, behavinsight_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, narmism_avg); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))


%%
% Create line plot with channel x timewindow trial data
meanLine = std(sl_data_test.trial,1);

figure('Renderer', 'painters', 'Position', [10 10 500 200])
hold on
plot(sl_data_test.time, sl_data_test.trial, 'color', [0.5 0.5 0.5]);
plot(sl_data_test.time, meanLine, 'color', [0.6350 0.0780 0.1840]);
ylim([-0.012 0.012])
xlim([1 15])
xlabel('Time-window')
ylabel('Similarity (r)')
title('Similarity for each channel over time')
set(gcf,'color','w')


% Create image of prediction matrix

predictionRDM(isnan(predictionRDM)) = 0;
imagesc(predictionRDM);
xlabel('Video Event')
ylabel('Video Event')
title('Prediction RSM - Behavior')
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 


%% Test out stuff

test = zeros(12195,15);


load(fullfile('/project/3012026.13/processed_RT/', 'searchlight_test.mat'))


%% Plot results

subj = subjects{1};

load(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'dataRSA.mat'));
dataRSAv1 = dataRSA; % comb planar, no filter
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'dataRSAv2.mat'));
dataRSAv2 = dataRSA; % axial grad, 35Hz filter
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'dataRSAv3.mat'));
dataRSAv3 = dataRSA; % axial grad, no filter
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'dataRSA_planarsplit.mat'));
dataRSAv4 = dataRSA; % planar split, no filter
load(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'dataRSA_planarsplitsubtrmean.mat'));
dataRSAv5 = dataRSA; % planar split, filter, mean subtracted

load(fullfile('/project/3012026.13/processed_RT/', 'template_timelock.mat'));
cfg                 = [];
cfg.avgoverrpt      = 'yes';
cfg.avgovertime     = 'no';
template_data       = ft_selectdata(cfg, template_timelock);
template_data.time  = [1:15]*0.125;
template_data       = rmfield(template_data, 'trial');
template_data.cfg   = [];
clear template_timelock

dataRSAv1.time = template_data.time;
dataRSAv1.elec = template_data.elec;
dataRSAv1.grad = template_data.grad;
dataRSAv1.dimord = template_data.dimord;
dataRSAv1 = rmfield(dataRSAv1, 'windownr_time');

dataRSAv2.time = template_data.time;
dataRSAv2.elec = template_data.elec;
dataRSAv2.grad = template_data.grad;
dataRSAv2.dimord = template_data.dimord;
dataRSAv2 = rmfield(dataRSAv2, 'windownr_time');

dataRSAv3.time = template_data.time;
dataRSAv3.elec = template_data.elec;
dataRSAv3.grad = template_data.grad;
dataRSAv3.dimord = template_data.dimord;
dataRSAv3 = rmfield(dataRSAv3, 'windownr_time');

dataRSAv4.time = template_data.time;
dataRSAv4.elec = template_data.elec;
dataRSAv4.grad = template_data.grad;
dataRSAv4.dimord = template_data.dimord;
dataRSAv4.label = template_data.label;
dataRSAv4 = rmfield(dataRSAv4, 'windownr_time');

dataRSAv5.time = template_data.time;
dataRSAv5.elec = template_data.elec;
dataRSAv5.grad = template_data.grad;
dataRSAv5.dimord = template_data.dimord;
dataRSAv5.label = template_data.label;
dataRSAv5 = rmfield(dataRSAv5, 'windownr_time');

% Topographies
% Plot the data onto topography
cfg                 = [];
cfg.xlim            = [1:15]*0.125;
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'vis1';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, dataRSAv1); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, dataRSAv2); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, dataRSAv3); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, dataRSAv4); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

figure; ft_topoplotER(cfg, dataRSAv5); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

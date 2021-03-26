%% Apply RSA with spatial and temporal searchlight to timelock data
%% Default path

addpath /home/common/matlab/fieldtrip/
addpath /home/common/matlab/fieldtrip/qsub/
addpath /project/3012026.13/scripts_RT/
addpath /project/3012026.13/scripts_RT/Scripts_MCCA_Module
ft_defaults

root_dir     = '/project/3012026.13/processed/';
save_dir     = '/project/3012026.13/processed_RT/time_lock/';

subjects = {'sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-022','sub-023',...
    'sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-031'...
    'sub-032','sub-033','sub-034','sub-035','sub-036','sub-037','sub-038'};


%% Load in timelock data and apply spatial/temporal searchlight

% These subjects are ignored because of inconsistend trial structures
ignSub = [4, 10, 12, 20, 30, 32];

% redo sub: 12(9) 20(16)

RSA_error = [];
get_neighbours = true;
for iSubject = 1:length(subjects)
    
    if ~ismember(iSubject, ignSub)

        disp(strcat('Initiating searchligh subject:', int2str(iSubject)))

        % Load Data
        subj                = subjects{iSubject};
        % Apply timelock analysis
        my_data             = rt_mytimelockv2(root_dir, subj);

        % Get neighbour structure, do it only once
        if get_neighbours
            load(fullfile('/project/3012026.13/processed_RT/', 'template_timelock.mat'));
            cfg             = [];
            cfg.method      = 'triangulation';
            cfg.layout      = 'CTF275_helmet.mat';
            cfg.channel     = template_timelock.label;
            neighbours      = ft_prepare_neighbours(cfg);
            clear template_timelock
        end
        
        try
            % Apply searchlight function
            cfg                 = [];
            cfg.timewidth       = 0.250; % Width of shifting timewindow
            cfg.timesteps       = 0.125; % Steps of shifting timewindow
            cfg.neighbours      = neighbours;
            cfg.searchspace     = 'yes'; % Not included in rt_searchlight yet
            cfg.searchtime      = 'yes'; % Not included in rt_searchlight yet
            cfg.avgovertrials   = 'yes'; % Not included in rt_searchlight yet
            cfg.avgovertime     = 'no'; % Not included in rt_searchlight yet
            cfg.avgoverspace    = 'no'; % Not included in rt_searchlight yet
            
            qsubfeval('rt_searchlight', cfg, my_data, subj, 'memreq', (1024^3)*12, 'timreq', 60*180);

            %rt_searchlight(cfg, my_data, subj);
        catch
            disp(strcat('Error sub:', int2str(iSubject)))
            RSA_error = [RSA_error, iSubject];
        end

        get_neighbours = false;

        keep iSubject subjects root_dir save_dir neighbours get_neighbours ignSub RSA_error

    end
    
end


%% Create prediction RDM and compare to neural RDM

%predata     = xlsread('/project/3012026.13/scripts_RT/prediction_pre_noBX.xlsx');
%postdata    = xlsread('/project/3012026.13/scripts_RT/prediction_post_noBX.xlsx');
%sanitydata  = xlsread('/project/3012026.13/scripts_RT/sanity_pre.xlsx');
%Create prediction RDM that matches trial x trial structure of data
cfg             = [];
cfg.offdiag     = false; % Use postdata structure for off diagional (phase comparison)
predictionRDM   = rt_predictionRDMxlsx(cfg, predata, postdata);

ignoreSubs = [4 10 12 20];
for subby = 1:20
    
  
    subj            = subjects{subby};
    subj            = strcat(subj(1:3), subj(5:length(subj)));
    sub_behav_path  = fullfile('/project/3012026.13/processed_RT/logfiles/', [subj,'_Relatedness_answers.txt']);
    predictionRDM   = rt_predictionRDMbehav(sub_behav_path);
    
    disp(strcat('Loading Subject:', int2str(subby)))
    load(fullfile('/project/3012026.13/processed_RT/searchlight/', strcat(subjects{subby},'searchlight.mat')))
    neuralRDM       = sl_data.searchlight;
    clear sl_data

    predictionRDM_vec   = vectorizeSimmat(predictionRDM);
    dataRSA{subby,1}             = zeros(size(neuralRDM,1), size(neuralRDM,2));
    
    if ~ismember(subby, ignoreSubs)
        for nrChan = 1:size(neuralRDM, 1)
            disp(strcat('Calculating Channel:', int2str(nrChan)))
            for nrTimewin = 1:size(neuralRDM, 2)

                currentNeural               = squeeze(neuralRDM(nrChan, nrTimewin, :, :));
                currentNeural               = vectorizeSimmat(currentNeural);
                % dataRSA is the resulting structure containing chan x timewindow
                % information
                dataRSA{subby,1}(nrChan, nrTimewin)  = corr(predictionRDM_vec', currentNeural', 'Type', 'Kendall', 'Rows', 'complete');

            end
        end
    end
    clear neuralRDM
end

% Throw out sub 4, 10, 12, 20
% First save original data
dataRSAoriginal = dataRSA;

dataRSA{4} = [];
dataRSA{10} = [];
dataRSA{12} = [];
dataRSA{20} = [];
dataRSA = dataRSA(~cellfun('isempty',dataRSA));

% Calculate the mean correlation over subjects
sumRSA = zeros(size(dataRSA{1},1), size(dataRSA{1},2));
for meanSub = 1:length(dataRSA)  
    sumRSA = sumRSA + dataRSA{meanSub};
end
meanRSA = sumRSA./length(dataRSA);

% Integrate data from RSA into structure that is intepretable by fieldtrip
load(fullfile('/project/3012026.13/processed_RT/', 'template_timelock.mat'));
cfg                 = [];
cfg.avgoverrpt      = 'yes';
cfg.avgovertime     = 'no';
sl_data_test        = ft_selectdata(cfg, template_timelock);
sl_data_test.trial  = meanRSA(:,:);
sl_data_test.time   = [1:15];

% Plot the data onto topography
cfg                 = [];
cfg.zlim            = [-0.012 0.012]; 
cfg.xlim            = [1:1:15];
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'trial';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
figure; ft_topoplotER(cfg,sl_data_test); colorbar
set(gcf,'color','w')

ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu'))) 

% Plot channel x timewindow data
figure('Renderer', 'painters', 'Position', [10 10 300 500])
imagesc(meanRSA, [-0.012 0.012])
xlabel('Timewindow')
ylabel('Channel')
title('Similarity (r)')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 


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

%% RSA based on MCCA

% Run for each subject
for m = 1:numel(subjects)
    subj = subjects{m};
    %qsubfeval('rt_mcca', subj, 'hilbert', 'abs', 'memreq', (1024^3)*12, 'timreq', 60*60);
    %qsubfeval('rt_mcca', subj, 'memreq', (1024^3)*12, 'timreq', 60*60); % rt_mytimelockv3 includes hilbert
    %qsubfeval('rt_mcca_ERP', subj, 'memreq', (1024^3)*12, 'timreq', 60*60); 
    qsubfeval('rt_mcca_FREQ', subj, 'memreq', (1024^3)*12, 'timreq', 60*60); 
end

% Collect data for group avg
%suff = 'mcca_theta_rsm';
%suff = 'mcca_alpha_rsm';
%suff = 'mcca_ERP';
%suff = 'mcca_FREQ';
suff = 'tfr.mat';
%rt_collect_rsm(suff)
%rt_collect_FREQ(suff)
rt_collectdata(suff, 'theta')
rt_collectdata(suff, 'alpha')
rt_collectdata(suff, 'beta')

% Load the aggregated data
load('/project/3012026.13/jansch/groupdata_mcca_theta_hilbert.mat')
load('/project/3012026.13/jansch/groupdata_mcca_alpha_hilbert.mat')
load('/project/3012026.13/jansch/groupdata_mcca_alpha_rsm.mat')
load('/project/3012026.13/jansch/groupdata_mcca_theta_rsm.mat')
load('/project/3012026.13/jansch/groupdata_mcca_ERP.mat')
load('/project/3012026.13/jansch/groupdata_mcca_FREQ.mat')


%% Do correlation of neural RSM and prediction RSM

suff = '_mcca_theta_rsm';
%datadir = '/project/3012026.13/jansch/subject_jobs';
%cd(datadir);
for m = 1:numel(subjects)
    m
    subj = subjects{m};
    %qsubfeval('rt_corr_mcca', subj, suff, 'memreq', (1024^3)*12, 'timreq', 60*60);
    rt_corr_mcca(subj, suff)
    disp('Artificial pause...')
    pause(720)
end

% Combine the individual channel data
for m = 1:numel(subjects)
    try 
        m
        subj = subjects{m};
        suff = 'mcca_chan';
        rt_combine_rsmChan(subj,suff)
    end
end


%% Load combined data and create new structure from template

load('/project/3012026.13/jansch/sub-004_mcca_theta_rsm.mat')
templatedata = tlck{1,1};


for m = 1:numel(subjects)
    m
    subj = subjects{m};
    try
        load(fullfile('/project/3012026.13/jansch/', subj, ['combChan_', subj, '.mat']))
        tlck = templatedata;
        tlck.avg = combChan;
        save(fullfile('/project/3012026.13/jansch/', [subj, '_predRSM', 'theta.mat']), 'tlck');
        clear tlck combChan
    end
end

% Get avg over subjects
datadir = '/project/3012026.13/jansch';
cd(datadir);
d = dir(sprintf('sub*%s*','predRSMtheta'));
for k = 1:numel(d)
  k
  load(d(k).name);
  groupPredRSM{k,1} = tlck;
end

cfg = [];
avgPredRSM = ft_timelockgrandaverage(cfg, groupPredRSM{:,1});


%% Plot

load('/project/3012026.13/jansch/groupdata_tlck.mat')
tlck_stat = stat;
load('/project/3012026.13/jansch/groupdata_theta_tfr.mat')
theta_stat = stat;
load('/project/3012026.13/jansch/groupdata_alpha_tfr.mat')
alpha_stat = stat;
clear stat

cfg = [];
cfg.layout = 'CTF275_helmet.mat';  
cfg.parameter = 'stat';
cfg.xlim = [0.3 0.43];
%ft_topoplotER(cfg, tlck_stat{1});
ft_topoplotTFR(cfg, theta_stat{1});

%% Plot RSA based on MCCA results

%    1 {'Apre'}    {'Apre'}
%    2 {'Apre'}    {'Bpre'}
%    3 {'Apre'}    {'Xpre'}
%    4 {'Apre'}    {'Apst'}
%    5 {'Apre'}    {'Bpst'}
%    6 {'Apre'}    {'Xpst'}
%    7 {'Bpre'}    {'Bpre'}
%    8 {'Bpre'}    {'Xpre'}
%    9 {'Bpre'}    {'Apst'}
%   10 {'Bpre'}    {'Bpst'}
%   11 {'Bpre'}    {'Xpst'}
%   12 {'Xpre'}    {'Xpre'}
%   13 {'Xpre'}    {'Apst'}
%   14 {'Xpre'}    {'Bpst'}
%   15 {'Xpre'}    {'Xpst'}
%   16 {'Apst'}    {'Apst'}
%   17 {'Apst'}    {'Bpst'}
%   18 {'Apst'}    {'Xpst'}
%   19 {'Bpst'}    {'Bpst'}
%   20 {'Bpst'}    {'Xpst'}
%   21 {'Xpst'}    {'Xpst'}

% Index for cmb matrix to transform into 6x6 matrix
matCMB = [1, 2, 3, 4, 5, 6; 2, 7, 8, 9, 10, 11; 3, 8, 12, 13, 14, 15;
    4, 5, 6, 16, 17, 18; 9, 10, 11, 17, 19, 20; 13, 14, 15, 18, 20, 21];

% ft_timelockgrandaverage needs avg field
for k = 1:numel(T)
  T{k}.avg = T{k}.trial;
end

MCCAavg = cell(size(T,1),1);
cfg = [];
for x = 1:size(T,1)
    MCCAavg{x,1} = ft_timelockgrandaverage(cfg, T{x,:});
end

% Devide data in time-windows similar to previous attempt
fsample = round(1./mean(diff(MCCAavg{1}.time))); % assume integer
nwidth = round(0.25.*fsample);
nstep = round(0.125.*fsample);
nmax = 15; % hardcoded to match previous attempt, i.e. 15 time-windows

% Create array that contains cellcmb X channel X timewindow data
MCCAtw = zeros(size(MCCAavg, 1), size(MCCAavg{1}.label, 1), nmax);
for s = 1:length(MCCAavg)
    twidth = nwidth;
    tstep = 1;
    for t = 1:nmax
        MCCAtw(s, :, t) = mean(MCCAavg{s}.avg(:, tstep:twidth), 2);
        tstep = tstep + nstep;
        twidth = twidth + nstep;
    end
end

% Create array that contains channel X time-window X 6x6 matrix
MCCAmat = zeros(size(MCCAavg{1}.label, 1), nmax, size(matCMB,1), size(matCMB,2));
for t = 1:size(MCCAtw, 3)
    for mr = 1:size(matCMB, 1) % for rows
        for mc = 1:size(matCMB, 2) % for columns
            MCCAmat(:, t, mr, mc) = squeeze(MCCAtw(matCMB(mr, mc), :, t));
        end
    end
end


%% Plot results

% Devide the brain into frontal,L/R temporal and posterior channels
% Plot RSMs per time window

antchan = {'MLC12', 'MLC13', 'MLC14', 'MLC15', 'MLC21', 'MLC22', 'MLC23', 'MLC24', 'MLC31', 'MLC41', 'MLC51', 'MLC52', 'MLC53', 'MLC61', 'MLC62', 'MLF23', 'MLF31', 'MLF32', 'MLF41', 'MLF42', 'MLF43', 'MLF51', 'MLF52', 'MLF53', 'MLF61', 'MLF62', 'MLF63', 'MRC11', 'MRC12', 'MRC13', 'MRC21', 'MRC22', 'MRC23', 'MRC31', 'MRC41', 'MRC51', 'MRC52', 'MRC53', 'MRC61', 'MRC62', 'MRF31', 'MRF32', 'MRF41', 'MRF42', 'MRF43', 'MRF51', 'MRF52', 'MRF61', 'MRF62', 'MZC01', 'MZC02', 'MZC03', 'MZF02', 'MZF03'};
postchan = {'MLO11', 'MLO12', 'MLO13', 'MLO21', 'MLO22', 'MLO23', 'MLO24', 'MLO31', 'MLO32', 'MLO41', 'MLO42', 'MLO43', 'MLP21', 'MLP31', 'MLP32', 'MLP33', 'MLP41', 'MLP42', 'MLP43', 'MLP51', 'MLP52', 'MLP53', 'MLP54', 'MRO11', 'MRO12', 'MRO13', 'MRO14', 'MRO21', 'MRO22', 'MRO23', 'MRO24', 'MRO31', 'MRO32', 'MRO33', 'MRO41', 'MRO42', 'MRO43', 'MRO44', 'MRO53', 'MRP21', 'MRP31', 'MRP32', 'MRP33', 'MRP41', 'MRP42', 'MRP43', 'MRP44', 'MRP51', 'MRP52', 'MRP53', 'MRP54', 'MRP55', 'MZO01', 'MZO02', 'MZP01'};
leftchan = {'MLC15', 'MLC16', 'MLC17', 'MLC24', 'MLC25', 'MLF35', 'MLF44', 'MLF45', 'MLF46', 'MLF54', 'MLF55', 'MLF56', 'MLF63', 'MLF64', 'MLF65', 'MLF66', 'MLF67', 'MLO14', 'MLP35', 'MLP43', 'MLP44', 'MLP45', 'MLP54', 'MLP55', 'MLP56', 'MLP57', 'MLT11', 'MLT12', 'MLT13', 'MLT14', 'MLT15', 'MLT16', 'MLT21', 'MLT22', 'MLT23', 'MLT24', 'MLT25', 'MLT26', 'MLT31', 'MLT32', 'MLT33', 'MLT34', 'MLT35', 'MLT36', 'MLT41', 'MLT42', 'MLT43', 'MLT44', 'MLT45', 'MLT53'};
rightchan = {'MRC13', 'MRC14', 'MRC15', 'MRC16', 'MRC17', 'MRC22', 'MRC23', 'MRC24', 'MRC25', 'MRC31', 'MRC32', 'MRC42', 'MRF25', 'MRF34', 'MRF35', 'MRF44', 'MRF45', 'MRF46', 'MRF52', 'MRF53', 'MRF54', 'MRF55', 'MRF56', 'MRF62', 'MRF63', 'MRF64', 'MRF65', 'MRF67', 'MRO14', 'MRP23', 'MRP33', 'MRP34', 'MRP35', 'MRP42', 'MRP43', 'MRP44', 'MRP45', 'MRP53', 'MRP54', 'MRP55', 'MRP56', 'MRP57', 'MRT11', 'MRT12', 'MRT13', 'MRT14', 'MRT15', 'MRT16', 'MRT21', 'MRT22', 'MRT23', 'MRT24', 'MRT25', 'MRT26', 'MRT27', 'MRT31', 'MRT32', 'MRT33', 'MRT34', 'MRT35', 'MRT36', 'MRT42', 'MRT43', 'MRT44', 'MRT45', 'MRT53', 'MRT54'};

% Anterior channels
mychannels  = antchan;
channels    = cell2mat(MCCAavg{1}.label);
dataselect  = MCCAmat(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    caxis([-0.25 0.25]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('MCCA: Anterior Electrodes Neural RSMs')

% Posterior channels
mychannels  = postchan;
channels    = cell2mat(MCCAavg{1}.label);
dataselect  = MCCAmat(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    caxis([-0.25 0.25]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('MCCA: Posterior Electrodes Neural RSMs')

% left channels
mychannels  = leftchan;
channels    = cell2mat(MCCAavg{1}.label);
dataselect  = MCCAmat(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    caxis([-0.25 0.25]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('MCCA: Left Lateral Electrodes Neural RSMs')

% left channels
mychannels  = rightchan;
channels    = cell2mat(MCCAavg{1}.label);
dataselect  = MCCAmat(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    caxis([-0.25 0.25]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('MCCA: Right Lateral Electrodes Neural RSMs')


%% Plot the difference between cells (pre-post)

% Calculate difference

neuralSim_post_pre_diff = MCCAmat(:,:,4:6,4:6) - MCCAmat(:,:,1:3,1:3);

% Anterior channels
mychannels  = antchan;
channels    = cell2mat(MCCAavg{1}.label);
dataselect  = neuralSim_post_pre_diff(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

% Barplot
figure
for images = 1:size(dataselect,1)
    ax(images)  = subplot(5,3,images);
    withinEvent = diag(squeeze(dataselect(images,:,:)))'; % AA, BB, XX
    acrossEvent = vectorizeSimmat(squeeze(dataselect(images,:,:))); % AB, AX, BX
    allComp     = [withinEvent, acrossEvent];
    bar(allComp)
    set(gca, 'XTickLabel', {'AA' 'BB' 'XX' 'AB' 'AX' 'BX'})
    set(gcf, 'Position',  [100, 100, 500, 800])
    xtickangle(45)
    ylim([-9*10^(-3) 9*10^(-3)])
    title(strcat('T', num2str(images)))
end
sgtitle('MCCA: Anterior Electrodes Post-Pre Phase Similarity')

% Posterior channels
mychannels  = postchan;
channels    = cell2mat(MCCAavg{1}.label);
dataselect  = neuralSim_post_pre_diff(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

% Barplot
figure
for images = 1:size(dataselect,1)
    ax(images)  = subplot(5,3,images);
    withinEvent = diag(squeeze(dataselect(images,:,:)))'; % AA, BB, XX
    acrossEvent = vectorizeSimmat(squeeze(dataselect(images,:,:))); % AB, AX, BX
    allComp     = [withinEvent, acrossEvent];
    bar(allComp)
    set(gca, 'XTickLabel', {'AA' 'BB' 'XX' 'AB' 'AX' 'BX'})
    set(gcf, 'Position',  [100, 100, 500, 800])
    xtickangle(45)
    ylim([-9*10^(-3) 9*10^(-3)])
    title(strcat('T', num2str(images)))
end
sgtitle('MCCA: Posterior Electrodes Post-Pre Phase Similarity')

% Left channels
mychannels  = leftchan;
channels    = cell2mat(MCCAavg{1}.label);
dataselect  = neuralSim_post_pre_diff(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

% Barplot
figure
for images = 1:size(dataselect,1)
    ax(images)  = subplot(5,3,images);
    withinEvent = diag(squeeze(dataselect(images,:,:)))'; % AA, BB, XX
    acrossEvent = vectorizeSimmat(squeeze(dataselect(images,:,:))); % AB, AX, BX
    allComp     = [withinEvent, acrossEvent];
    bar(allComp)
    set(gca, 'XTickLabel', {'AA' 'BB' 'XX' 'AB' 'AX' 'BX'})
    set(gcf, 'Position',  [100, 100, 500, 800])
    xtickangle(45)
    ylim([-9*10^(-3) 9*10^(-3)])
    title(strcat('T', num2str(images)))
end
sgtitle('MCCA: Left Lateral Electrodes Post-Pre Phase Similarity')

% Right channels
mychannels  = rightchan;
channels    = cell2mat(MCCAavg{1}.label);
dataselect  = neuralSim_post_pre_diff(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

% Barplot
figure
for images = 1:size(dataselect,1)
    ax(images)  = subplot(5,3,images);
    withinEvent = diag(squeeze(dataselect(images,:,:)))'; % AA, BB, XX
    acrossEvent = vectorizeSimmat(squeeze(dataselect(images,:,:))); % AB, AX, BX
    allComp     = [withinEvent, acrossEvent];
    bar(allComp)
    set(gca, 'XTickLabel', {'AA' 'BB' 'XX' 'AB' 'AX' 'BX'})
    set(gcf, 'Position',  [100, 100, 500, 800])
    xtickangle(45)
    ylim([-9*10^(-3) 9*10^(-3)])
    title(strcat('T', num2str(images)))
end
sgtitle('MCCA: Right Lateral Electrodes Post-Pre Phase Similarity')


%% Plot from stats channels

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'CTF275_helmet.mat';
cfg.channel = stat_AXpost_AXpre.label(idxProb);
ft_singleplotER(cfg, MCCAavg{18}, MCCAavg{3});



%%
post = MCCAavg{18};
pre = MCCAavg{3};
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
ft_multiplotER(cfg, post, pre)
















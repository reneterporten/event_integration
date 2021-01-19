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


%% Collect all data across subjects


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
        
        % Create data structure from template ERP data
        dataRSApre                  = template_data;
        dataRSApost                 = template_data;
        
        dataRSApre.narinsight   = dataRSA.narinsight_pre;
        dataRSApost.narinsight   = dataRSA.narinsight_post;
        
        % Visual similarity
%         dataRSAorg.vis1             = dataRSA.vis1;
%         dataRSAperm.vis1            = dataRSA.vis1_rand;
%         
%         % Narrative insight interaction
%         dataRSAorg.narinsight       = dataRSA.narinsight;
%         dataRSAperm.narinsight      = dataRSA.narinsight_rand;
%         
%         % Narrative insight behavior
%         dataRSAorg.behavinsight     = dataRSA.behavinsight;
%         dataRSAperm.behavinsight    = dataRSA.behavinsight_rand;
%         
%         % Narrative mismatch
%         dataRSAorg.narmism          = dataRSA.narmism;
%         dataRSAperm.narmism         = dataRSA.narmism_rand;
%         
%         % (Un)lining mechanism
%         dataRSAorg.unlink           = dataRSA.unlink;
%         dataRSAperm.unlink          = dataRSA.unlink_rand;

        allDataRSApre{rsa_count}    = dataRSApre;
        allDataRSApost{rsa_count}   = dataRSApost;
        rsa_count                   = rsa_count + 1;
    
    end

end

keep allDataRSApost allDataRSApre ignSub subjects save_dir root_dir


%% Cluster based stats across subjects (within groups)

cfg                         = [];
cfg.method                  = 'template'; 
cfg.template                = 'CTF275_neighb.mat';               
cfg.layout                  = 'CTF275_helmet.mat';                     
cfg.feedback                = 'no';                            
neighbours                  = ft_prepare_neighbours(cfg, allDataRSApre{1}); 

cfg                         = [];
cfg.channel                 = {'all'};
cfg.neighbours              = neighbours;
cfg.latency                 = 'all';
cfg.method                  = 'montecarlo';
cfg.statistic               = 'ft_statfun_depsamplesT';
cfg.correctm                = 'cluster';
cfg.avgovertime             = 'no';
cfg.clusteralpha            = 0.05;
cfg.clusterstatistic        = 'maxsum';
cfg.minnbchan               = 2;
cfg.tail                    = 0;
cfg.clustertail             = 0;
cfg.alpha                   = 0.05;
cfg.correcttail             = 'prob';
cfg.numrandomization        = 1000;

Nsubj  = length(allDataRSApre);
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

cfg.parameter       = 'narinsight';
stat_nar_prepost    = ft_timelockstatistics(cfg, allDataRSApre{:}, allDataRSApost{:});
% cfg.parameter       = 'vis1';
% stat_vis            = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});
% cfg.parameter       = 'narinsight';
% stat_narinsight     = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});
% cfg.parameter       = 'behavinsight';
% stat_behavinsight   = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});
% cfg.parameter       = 'narmism';
% stat_narmism        = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});
% cfg.parameter       = 'unlink';
% stat_unlink         = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});

%save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'stat_preeff_vs_posteff.mat'), 'stat_exp_vs_surp', '-v7.3')
%disp('Done saving.')    


%% Cluster based stats across subjects (between groups)

cfg                         = [];
cfg.method                  = 'template'; 
cfg.template                = 'CTF275_neighb.mat';               
cfg.layout                  = 'CTF275_helmet.mat';                     
cfg.feedback                = 'no';                            
neighbours                  = ft_prepare_neighbours(cfg, allDataRSAorg{1}); 

cfg                         = [];
cfg.channel                 = {'all'};
cfg.neighbours              = neighbours;
cfg.latency                 = 'all';
cfg.method                  = 'montecarlo';
cfg.statistic               = 'ft_statfun_indepsamplesT';
cfg.correctm                = 'cluster';
cfg.avgovertime             = 'no';
cfg.clusteralpha            = 0.05;
cfg.clusterstatistic        = 'maxsum';
cfg.minnbchan               = 2;
cfg.tail                    = 0;
cfg.clustertail             = 0;
cfg.alpha                   = 0.05;
cfg.correcttail             = 'prob';
cfg.numrandomization        = 1000;

subj = length(allDataRSAorg);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end 
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;


cfg.design = design;
cfg.ivar  = 2;

cfg.parameter       = 'vis1';
stat_vis            = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});
cfg.parameter       = 'narinsight';
stat_narinsight     = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});
cfg.parameter       = 'behavinsight';
stat_behavinsight   = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});
cfg.parameter       = 'narmism';
stat_narmism        = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});
cfg.parameter       = 'unlink';
stat_unlink         = ft_timelockstatistics(cfg, allDataRSAorg{:}, allDataRSAperm{:});

%save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'stat_preeff_vs_posteff.mat'), 'stat_exp_vs_surp', '-v7.3')
%disp('Done saving.')    


%% Stats on MCCA data

cfg                         = [];
cfg.method                  = 'template'; 
cfg.template                = 'CTF275_neighb.mat';               
cfg.layout                  = 'CTF275_helmet.mat';                     
cfg.feedback                = 'no';                            
neighbours                  = ft_prepare_neighbours(cfg, T{1,1}); 

cfg                         = [];
cfg.channel                 = {'all'};
cfg.neighbours              = neighbours;
cfg.latency                 = 'all';
cfg.method                  = 'montecarlo';
cfg.statistic               = 'ft_statfun_depsamplesT';
cfg.correctm                = 'cluster';
cfg.avgovertime             = 'no';
cfg.clusteralpha            = 0.05;
cfg.clusterstatistic        = 'maxsum';
cfg.minnbchan               = 2;
cfg.tail                    = 0;
cfg.clustertail             = 0;
cfg.alpha                   = 0.05;
cfg.correcttail             = 'prob';
cfg.numrandomization        = 1000;

Nsubj  = size(T,2);
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

stat_ABpost_ABpre    = ft_timelockstatistics(cfg, T{17,:}, T{2,:});
stat_AXpost_AXpre    = ft_timelockstatistics(cfg, T{18,:}, T{3,:});


cfg = []; cfg.parameter = 'stat'; cfg.layout = 'CTF275_helmet.mat'; ft_multiplotER(cfg, stat_AXpost_AXpre)


%% Cluster plot of stats

stat_vis.cfg = [];
stat_narinsight.cfg = [];
stat_behavinsight.cfg = [];
stat_narmism.cfg = [];
stat_unlink.cfg = [];

cfg = [];
cfg.alpha = 0.06;
cfg.layout = 'CTF275_helmet.mat'; 
ft_clusterplot(cfg, stat_AXpost_AXpre)
ft_hastoolbox('brewermap', 1);
set(gcf,'color','w')
colormap(flipud(brewermap(64,'RdBu')))

% Define channels to be highlighted
tstart = -0.19972;
tstop = 0.33714;
[val,idxstart] = min(abs(stat_AXpost_AXpre.time-tstart));
[val,idxstop] = min(abs(stat_AXpost_AXpre.time-tstop));

% Index of channels that are below specific p-val
idxProb = stat_AXpost_AXpre.prob(:,idxstart:idxstop) < 0.06;
idxProb = sum(idxProb,2);
idxProb = idxProb > 0;

cfg = [];
cfg.xlim = [-0.19972 0.33714];
cfg.parameter = 'stat';
cfg.layout = 'CTF275_helmet.mat';
cfg.highlight = 'on';
cfg.highlightchannel = stat_AXpost_AXpre.label(idxProb);
cfg.highlightsymbol = '.';
cfg.highlightsize = 24;
ft_topoplotER(cfg, stat_AXpost_AXpre)
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

cfg = [];
cfg.parameter = 'stat';
cfg.layout = 'CTF275_helmet.mat';
cfg.channel = stat_AXpost_AXpre.label(idxProb);
ft_singleplotER(cfg, stat_AXpost_AXpre);


%% Plot of thresholded t-values

% Narrative insight pre vs post
% Select top 10% values
currentstat         = stat_nar_prepost;
rangetop            = (length(currentstat.stat(:))/100)*10;
currentstat_top5    = currentstat.stat;
[~,idx]             = sort(currentstat_top5(:), 'descend');
currentstat_toppos  = currentstat_top5(idx(1:round(rangetop)));
[~,idx]             = sort(currentstat_top5(:), 'ascend');
currentstat_topneg  = currentstat_top5(idx(1:round(rangetop)));
currentstat_masked  = currentstat;
topmask             = ismember(currentstat.stat, currentstat_toppos) + ismember(currentstat.stat, currentstat_topneg);
topmask             = topmask > 0;
currentstat_masked.stat(~topmask) = 0;

plotrange = [abs(mean(currentstat_topneg)), abs(mean(currentstat_toppos))];
plotrange = max(plotrange);

cfg                 = [];
cfg.xlim            = [1:15]*0.125;
cfg.zlim            = [-1*plotrange, plotrange];
%cfg.maskparameter   = topmask;
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'stat';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, currentstat_masked); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))

% Visual similarity
% Select top 10% values
currentstat         = stat_vis;
rangetop            = (length(currentstat.stat(:))/100)*10;
currentstat_top5    = currentstat.stat;
[~,idx]             = sort(currentstat_top5(:), 'descend');
currentstat_toppos  = currentstat_top5(idx(1:round(rangetop)));
[~,idx]             = sort(currentstat_top5(:), 'ascend');
currentstat_topneg  = currentstat_top5(idx(1:round(rangetop)));
currentstat_masked  = currentstat;
topmask             = ismember(currentstat.stat, currentstat_toppos) + ismember(currentstat.stat, currentstat_topneg);
topmask             = topmask > 0;
currentstat_masked.stat(~topmask) = 0;

plotrange = [abs(mean(currentstat_topneg)), abs(mean(currentstat_toppos))];
plotrange = max(plotrange);

cfg                 = [];
%cfg.xlim            = [1:15]*0.125;
cfg.zlim            = [-1*plotrange, plotrange];
%cfg.maskparameter   = topmask;
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'stat';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, currentstat_masked); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))


% Interaction
% Select top 10% values
currentstat         = stat_narinsight;
rangetop            = (length(currentstat.stat(:))/100)*10;
currentstat_top5    = currentstat.stat;
[~,idx]             = sort(currentstat_top5(:), 'descend');
currentstat_toppos  = currentstat_top5(idx(1:round(rangetop)));
[~,idx]             = sort(currentstat_top5(:), 'ascend');
currentstat_topneg  = currentstat_top5(idx(1:round(rangetop)));
currentstat_masked  = currentstat;
topmask             = ismember(currentstat.stat, currentstat_toppos) + ismember(currentstat.stat, currentstat_topneg);
topmask             = topmask > 0;
currentstat_masked.stat(~topmask) = 0;

plotrange = [abs(mean(currentstat_topneg)), abs(mean(currentstat_toppos))];
plotrange = mean(plotrange);

cfg                 = [];
%cfg.xlim            = [1:15]*0.125;
cfg.zlim            = [-1*plotrange, plotrange];
%cfg.maskparameter   = topmask;
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'stat';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, currentstat_masked); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))


% Behavior
% Select top 10% values
currentstat         = stat_behavinsight;
rangetop            = (length(currentstat.stat(:))/100)*10;
currentstat_top5    = currentstat.stat;
[~,idx]             = sort(currentstat_top5(:), 'descend');
currentstat_toppos  = currentstat_top5(idx(1:round(rangetop)));
[~,idx]             = sort(currentstat_top5(:), 'ascend');
currentstat_topneg  = currentstat_top5(idx(1:round(rangetop)));
currentstat_masked  = currentstat;
topmask             = ismember(currentstat.stat, currentstat_toppos) + ismember(currentstat.stat, currentstat_topneg);
topmask             = topmask > 0;
currentstat_masked.stat(~topmask) = 0;

plotrange = [abs(mean(currentstat_topneg)), abs(mean(currentstat_toppos))];
plotrange = mean(plotrange);

cfg                 = [];
%cfg.xlim            = [1:15]*0.125;
cfg.zlim            = [-1*plotrange, plotrange];
%cfg.zlim            = [-1.0, 1.0];
%cfg.maskparameter   = topmask;
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'stat';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, currentstat_masked); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))


% Narrative mismatch
% Select top 10% values
currentstat         = stat_narmism;
rangetop            = (length(currentstat.stat(:))/100)*10;
currentstat_top5    = currentstat.stat;
[~,idx]             = sort(currentstat_top5(:), 'descend');
currentstat_toppos  = currentstat_top5(idx(1:round(rangetop)));
[~,idx]             = sort(currentstat_top5(:), 'ascend');
currentstat_topneg  = currentstat_top5(idx(1:round(rangetop)));
currentstat_masked  = currentstat;
topmask             = ismember(currentstat.stat, currentstat_toppos) + ismember(currentstat.stat, currentstat_topneg);
topmask             = topmask > 0;
currentstat_masked.stat(~topmask) = 0;

plotrange = [abs(mean(currentstat_topneg)), abs(mean(currentstat_toppos))];
plotrange = mean(plotrange);

cfg                 = [];
%cfg.xlim            = [1:15]*0.125;
cfg.zlim            = [-1*plotrange, plotrange];
%cfg.zlim            = [-1.0, 1.0];
%cfg.maskparameter   = topmask;
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'stat';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, currentstat_masked); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))


% (Un)linking mechanism
% Select top 10% values
currentstat         = stat_unlink;
rangetop            = (length(currentstat.stat(:))/100)*10;
currentstat_top5    = currentstat.stat;
[~,idx]             = sort(currentstat_top5(:), 'descend');
currentstat_toppos  = currentstat_top5(idx(1:round(rangetop)));
[~,idx]             = sort(currentstat_top5(:), 'ascend');
currentstat_topneg  = currentstat_top5(idx(1:round(rangetop)));
currentstat_masked  = currentstat;
topmask             = ismember(currentstat.stat, currentstat_toppos) + ismember(currentstat.stat, currentstat_topneg);
topmask             = topmask > 0;
currentstat_masked.stat(~topmask) = 0;

plotrange = [abs(mean(currentstat_topneg)), abs(mean(currentstat_toppos))];
plotrange = mean(plotrange);

cfg                 = [];
%cfg.xlim            = [1:15]*0.125;
cfg.zlim            = [-1*plotrange, plotrange];
%cfg.zlim            = [-1.0, 1.0];
%cfg.maskparameter   = topmask;
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'stat';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, currentstat_masked); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))
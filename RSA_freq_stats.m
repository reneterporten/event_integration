%% Apply RSA with spatial and temporal searchlight to time frequency data
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



%% Prepare differences for interaction

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
Fdiff_preAX = cell(1,length(F));
Fdiff_preBX = cell(1,length(F));
for subs = 1:length(F)
    Fdiff_preAX{subs} = ft_math(cfg, F{1,subs}, F{3,subs});
    Fdiff_preBX{subs} = ft_math(cfg, F{2,subs}, F{3,subs});
end

Fdiff_postAX = cell(1,length(F));
Fdiff_postBX = cell(1,length(F));
for subs = 1:length(F)
    Fdiff_postAX{subs} = ft_math(cfg, F{4,subs}, F{6,subs});
    Fdiff_postBX{subs} = ft_math(cfg, F{5,subs}, F{6,subs});
end

cfg = [];
cfg.parameter = 'avg';
Fdiff_preAX_avg = ft_freqgrandaverage(cfg, Fdiff_preAX{:});

figure; cfg = []; cfg.xlim = [0.5 1.0]; cfg.ylim = [15 20]; cfg.parameter = 'avg'; cfg.layout = 'CTF275_helmet.mat'; ft_topoplotTFR(cfg, Fdiff_preAX_avg)


%% Stats on MCCA data

cfg                         = [];
cfg.method                  = 'template'; 
cfg.template                = 'CTF275_neighb.mat';               
cfg.layout                  = 'CTF275_helmet.mat';                     
cfg.feedback                = 'no';                            
neighbours                  = ft_prepare_neighbours(cfg, F{1,1}); 

cfg                         = [];
cfg.channel                 = {'all'};
cfg.parameter               = 'avg';
cfg.neighbours              = neighbours;
cfg.frequency               = [8 12];
cfg.avgoverfreq             = 'yes';
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

Nsubj  = size(F,2);
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

% Testing within a 2x2 factor within subject anova framework
stat_interactionAX    = ft_freqstatistics(cfg, Fdiff_preAX{:}, Fdiff_postAX{:});
stat_interactionBX    = ft_freqstatistics(cfg, Fdiff_preBX{:}, Fdiff_postBX{:});


cfg = []; cfg.parameter = 'stat'; cfg.layout = 'CTF275_helmet.mat'; ft_multiplotER(cfg, stat_interactionAX)
cfg = []; cfg.parameter = 'stat'; cfg.alpha = 0.05; cfg.layout = 'CTF275_helmet.mat'; ft_clusterplot(cfg, stat_interactionAX)

cfg = []; cfg.xlim = [0.5 1.0]; cfg.parameter = 'stat'; cfg.layout = 'CTF275_helmet.mat'; ft_topoplotTFR(cfg, stat_interactionAX)

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
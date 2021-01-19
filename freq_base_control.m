%% Compare power spectra across control and baseline conditions

addpath /home/common/matlab/fieldtrip/
addpath /home/common/matlab/fieldtrip/qsub/
addpath /project/3012026.13/scripts_RT/
addpath /project/3012026.13/MVPA/MVPA-Light/startup
ft_defaults
startup_MVPA_Light

root_dir     = '/project/3012026.13/processed/';
save_dir     = '/project/3012026.13/processed_RT/time_freq/';

subjects = {'sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-022','sub-023',...
    'sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-031'...
    'sub-032','sub-033','sub-034','sub-035','sub-036','sub-037','sub-038'};


%% Time-freqeuncy queue

for iSubject = 1:length(subjects)

    subj        = subjects{iSubject};

    % load data, remove bad trials and demean & detrend
    load (fullfile(root_dir, [subj,'_dataclean.mat']));

    % remove bad trials
    load (fullfile(root_dir,'rejectedTrials', [subj,'_rejectedTrials.mat']));

    cfg         = [];
    cfg.trials  = not(rejected_logical);
    data        = ft_selectdata(cfg, dataclean);

    % detrend and demean the data
    cfg         = [];
    cfg.detrend = 'yes';
    data        = ft_preprocessing(cfg, data);

    % compute planar gradient
    cfg              = [];
    cfg.feedback     = 'no';
    cfg.method       = 'template';
    cfg.template     = 'CTF275_neighb.mat';
    cfg.planarmethod = 'sincos';
    cfg.channel      = {'MEG'};
    cfg.trials       = 'all';
    cfg.neighbours   = ft_prepare_neighbours(cfg, data);
    data_planar      = ft_megplanar(cfg,data);

    % select data
    % Baseline including
    A_pre_logical   = data_planar.trialinfo(:,2) == 1;
    A_post_logical  = data_planar.trialinfo(:,2) == 5;
    
    % Critical period
    B_pre_logical   = data.trialinfo(:,2) == 2;
    B_post_logical  = data.trialinfo(:,2) == 6;

    % Control
    X_pre_logical   = data_planar.trialinfo(:,2) == 3;
    X_post_logical  = data_planar.trialinfo(:,2) == 7;

    % Calculate tf spectra
    cfg             = [];
    cfg.trials      = A_pre_logical;
    data_A_pre      = ft_selectdata(cfg, data_planar);

    cfg             = [];
    cfg.trials      = A_post_logical;
    data_A_post     = ft_selectdata(cfg, data_planar);
    
    cfg.trials      = B_pre_logical;
    data_B_pre      = ft_selectdata(cfg, data_planar);

    cfg             = [];
    cfg.trials      = B_post_logical;
    data_B_post     = ft_selectdata(cfg, data_planar);

    cfg             = [];
    cfg.trials      = X_pre_logical;
    data_X_pre      = ft_selectdata(cfg, data_planar);

    cfg             = [];
    cfg.trials      = X_post_logical;
    data_X_post     = ft_selectdata(cfg, data_planar);

    qsubfeval('TF_basecont_qsubRT', save_dir, subj, data_A_pre, data_A_post, data_B_pre, data_B_post, data_X_pre, data_X_post, 'memreq', 24*(1024^3), 'timreq',1*60*60,'memoverhead', 24*(1024^3), 'batchid',sprintf('TF_%s',subj)); 

    keep iSubject subjects root_dir save_dir
    
end


%% Gather descriptives and calculate grand average contrast for X control

all_X_pre = {};
all_X_post = {};
for iSubject = 1:length(subjects)
    
    subj = subjects{iSubject};
    load (fullfile(save_dir, [subj,'_TF_basecon']), 'freq_X_pre','freq_X_post');
    
    cfg                     = [];
    cfg.keeptrials          = 'no';
    all_X_pre{iSubject}     = ft_freqdescriptives(cfg, freq_X_pre);
    all_X_post{iSubject}    = ft_freqdescriptives(cfg, freq_X_post);
    
    disp(int2str(iSubject))
    
    keep iSubject subjects root_dir save_dir all_X_pre all_X_post
     
end

freq_con_diff   = {};
for iSubject = 1:length(subjects)
    cfg                         = [];
    cfg.parameter               = 'powspctrm';
    cfg.operation               = 'subtract';
    freq_con_diff{iSubject}     = ft_math(cfg, all_X_pre{iSubject}, all_X_post{iSubject});
end

cfg             = [];
AVG_con_diff    = ft_freqgrandaverage(cfg, freq_con_diff{iSubject});


%% Gather descriptives and calculate grand average contrast for pre A baseline

all_A_pre = {};
all_A_post = {};
for iSubject = 1:length(subjects)
    
    subj = subjects{iSubject};
    load (fullfile(save_dir, [subj,'_TF_basecon']), 'freq_A_pre','freq_A_post');
    
    cfg                     = [];
    cfg.keeptrials          = 'no';
    all_A_pre{iSubject}     = ft_freqdescriptives(cfg, freq_A_pre);
    all_A_post{iSubject}    = ft_freqdescriptives(cfg, freq_A_post);
    
    disp(int2str(iSubject))
    
    keep iSubject subjects root_dir save_dir all_A_pre all_A_post
     
end

freq_base_diff   = {};
for iSubject = 1:length(subjects)
    cfg                         = [];
    cfg.parameter               = 'powspctrm';
    cfg.operation               = 'subtract';
    freq_base_diff{iSubject}     = ft_math(cfg, all_A_pre{iSubject}, all_A_post{iSubject});
end

cfg             = [];
AVG_base_diff    = ft_freqgrandaverage(cfg, freq_base_diff{iSubject});


%% Compare control conditions statistically 

cfg                     = [];
cfg.method              = 'template'; 
cfg.template            = 'CTF275_neighb.mat';               
cfg.layout              = 'CTF275_helmet.mat';                     
cfg.feedback            = 'no';                            
neighbours              = ft_prepare_neighbours(cfg, all_A_pre{1}); 

cfg                     = [];
cfg.neighbours          = neighbours;
cfg.channel             = 'all';

cfg.latency             = [-0.75 -0.25];
cfg.frequency           = [2 40];
cfg.avgovertime         = 'no';
cfg.avgoverfreq         = 'no';
cfg.avgoverchan         = 'no';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 2;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.025;
cfg.numrandomization    = 1000;

subj = length(all_A_pre);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end

design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;


cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;

%statX_pre_vs_post       = ft_freqstatistics(cfg,  all_A_pre{:},  all_A_post{:});
statX_pre_vs_post_time       = ft_freqstatistics(cfg,  all_A_pre{:},  all_A_post{:});

%%

cfg             = [];
cfg.alpha       = 0.05;
cfg.parameter   = 'stat';
cfg.layout      = 'CTF275_helmet.mat';
ft_clusterplot(cfg, statX_pre_vs_post)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu')))


%%
ind_sign_chan = find(statX_pre_vs_post.negclusterslabelmat == 1);
myChannels = {};
for ind=1:length(ind_sign_chan)
    
    myChannels{ind} = statX_pre_vs_post.label{ind_sign_chan(ind)};
    
end


cfg             = [];
cfg.channel     = myChannels;
%cfg.layout      = 'CTF275_helmet.mat';
cfg.xlim        = [-0.75 -0.25];
cfg.zlim        = [-2.5 2.5];
cfg.parameter   = 'stat';
figure
ft_singleplotTFR(cfg, statX_pre_vs_post_time);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu')))

%%

cfg             = [];
cfg.channel     = 'MLC12_dH';
cfg.layout      = 'CTF275_helmet.mat';
cfg.xlim        = [-0.75 0];
cfg.zlim        = [-3.0 3.0];
cfg.parameter   = 'stat';
figure
ft_multiplotTFR(cfg, statX_pre_vs_post );
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap


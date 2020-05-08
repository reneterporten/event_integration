%% Compare power spectra across control and baseline conditions

addpath /home/common/matlab/fieldtrip/
addpath /home/common/matlab/fieldtrip/qsub/
addpath /project/3012026.13/scripts_RT/
addpath /project/3012026.13/MVPA/MVPA-Light/startup
ft_defaults
startup_MVPA_Light

root_dir     = '/project/3012026.13/processed/';
save_dir     = '/project/3012026.13/processed_RT/time_lock/';

subjects = {'sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-022','sub-023',...
    'sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-031'...
    'sub-032','sub-033','sub-034','sub-035','sub-036','sub-037','sub-038'};


%% Test define trials

subj        = subjects{10};

% load data, remove bad trials and demean & detrend
load (fullfile(root_dir, [subj,'_dataclean.mat']));

cfg         = [];
cfg.trials  = find(dataclean.trialinfo(:,5)==1);
data_pre    = ft_selectdata(cfg, dataclean);

cfg         = [];
cfg.trials  = find(dataclean.trialinfo(:,5)==2);
data_post   = ft_selectdata(cfg, dataclean);

% Run through all the stories
data_pre_story  = {};
data_post_story = {};
for storyNr = 1:max(unique(data_post.trialinfo(:,4)))
    
    cfg         = [];
    cfg.trials  = find(data_pre.trialinfo(:,4) == storyNr);
    data_pre_story{storyNr,1} = ft_selectdata(cfg, data_pre);
    
    % It appears that the linking event has been coded in as well (first 6 trials of post)
    post_trials = find(data_post.trialinfo(:,4) == storyNr);
    post_trials = post_trials(7:end);
    
    cfg         = [];
    cfg.trials  = post_trials;
    data_post_story{storyNr,1} = ft_selectdata(cfg, data_post);
    
end


data_pre_story_A = {};
data_pre_story_B = {};
data_pre_story_X = {};

data_post_story_A = {};
data_post_story_B = {};
data_post_story_X = {};

for numStor = 1:size(data_pre_story,1)
    
    % Pre data
    cfg                         = [];
    cfg.trials                  = find(data_pre_story{numStor,1}.trialinfo(:,3) == 1);
    data_pre_story_A{numStor,1} = ft_selectdata(cfg, data_pre_story{numStor,1});
    
    cfg                         = [];
    cfg.trials                  = find(data_pre_story{numStor,1}.trialinfo(:,3) == 2);
    data_pre_story_B{numStor,1} = ft_selectdata(cfg, data_pre_story{numStor,1});
    
    cfg                         = [];
    cfg.trials                  = find(data_pre_story{numStor,1}.trialinfo(:,3) == 3);
    data_pre_story_X{numStor,1} = ft_selectdata(cfg, data_pre_story{numStor,1});
    
    % Post data
    cfg                             = [];
    cfg.trials                      = find(data_post_story{numStor,1}.trialinfo(:,3) == 1);
    data_post_story_A{numStor,1}    = ft_selectdata(cfg, data_post_story{numStor,1});
    
    cfg                             = [];
    cfg.trials                      = find(data_post_story{numStor,1}.trialinfo(:,3) == 2);
    data_post_story_B{numStor,1}    = ft_selectdata(cfg, data_post_story{numStor,1});
    
    cfg                             = [];
    cfg.trials                      = find(data_post_story{numStor,1}.trialinfo(:,3) == 3);
    data_post_story_X{numStor,1}    = ft_selectdata(cfg, data_post_story{numStor,1});
    
end


%% Time-lock analysis of the phase x story x video x trial data

subj        = subjects{10};

erp_story_pre   = {};
erp_story_post  = {};
for prepost = 1:size(data_pre_story_A,1)
    
    disp(strcat('Evaluating Story:', int2str(prepost)))
    
    % compute planar gradient
    cfg              = [];
    cfg.feedback     = 'no';
    cfg.method       = 'template';
    cfg.template     = 'CTF275_neighb.mat';
    cfg.planarmethod = 'sincos';
    cfg.channel      = {'MEG'};
    cfg.trials       = 'all';
    cfg.neighbours   = ft_prepare_neighbours(cfg, data_pre_story_A{prepost,1});
    
    data_planar_preA = ft_megplanar(cfg, data_pre_story_A{prepost,1});
    data_planar_preB = ft_megplanar(cfg, data_pre_story_B{prepost,1});
    data_planar_preX = ft_megplanar(cfg, data_pre_story_X{prepost,1});
    
    data_planar_postA = ft_megplanar(cfg, data_post_story_A{prepost,1});
    data_planar_postB = ft_megplanar(cfg, data_post_story_B{prepost,1});
    data_planar_postX = ft_megplanar(cfg, data_post_story_X{prepost,1});
    
    cfg             = [];
    cfg.channel     = 'MEG';
    cfg.trials      = 'all';
    cfg.keeptrials  = 'yes';
    erp_A_pre       = ft_timelockanalysis(cfg, data_planar_preA);
    erp_B_pre       = ft_timelockanalysis(cfg, data_planar_preB);
    erp_X_pre       = ft_timelockanalysis(cfg, data_planar_preX);
    erp_A_post      = ft_timelockanalysis(cfg, data_planar_postA);
    erp_B_post      = ft_timelockanalysis(cfg, data_planar_postB);
    erp_X_post      = ft_timelockanalysis(cfg, data_planar_postX);

    cfg             = [];
    cfg.method      = 'sum';
    erp_A_pre       = ft_combineplanar(cfg,erp_A_pre);
    erp_B_pre       = ft_combineplanar(cfg,erp_B_pre);
    erp_X_pre       = ft_combineplanar(cfg,erp_X_pre);
    erp_A_post      = ft_combineplanar(cfg,erp_A_post);
    erp_B_post      = ft_combineplanar(cfg,erp_B_post);
    erp_X_post      = ft_combineplanar(cfg,erp_X_post);

    cfg             = [];
    cfg.latency     = [0 2.0];
    erp_A_pre       = ft_selectdata(cfg,erp_A_pre);
    erp_B_pre       = ft_selectdata(cfg,erp_B_pre);
    erp_X_pre       = ft_selectdata(cfg,erp_X_pre);
    erp_A_post      = ft_selectdata(cfg,erp_A_post);
    erp_B_post      = ft_selectdata(cfg,erp_B_post);
    erp_X_post      = ft_selectdata(cfg,erp_X_post);
    
    erp_story_pre{prepost,1}.erp_A_pre   = erp_A_pre;
    erp_story_pre{prepost,1}.erp_B_pre   = erp_B_pre;
    erp_story_pre{prepost,1}.erp_X_pre   = erp_X_pre;
    
    erp_story_post{prepost,1}.erp_A_post   = erp_A_post;
    erp_story_post{prepost,1}.erp_B_post   = erp_B_post;
    erp_story_post{prepost,1}.erp_X_post   = erp_X_post;
    
end


%% Arrange data structure for RSA

my_data     = {};
storCounterA = 1;
storCounterB = 2;
storCounterX = 3;
for storEntries = 1:size(erp_story_pre,1)
    
    disp(strcat('Evaluating Story Pre:', int2str(storEntries)))
    
    my_data{storCounterA, 1} = erp_story_pre{storEntries,1}.erp_A_pre;
    my_data{storCounterB, 1} = erp_story_pre{storEntries,1}.erp_B_pre;
    my_data{storCounterX, 1} = erp_story_pre{storEntries,1}.erp_X_pre;
    
    storCounterA = storCounterA + 3;
    storCounterB = storCounterB + 3;
    storCounterX = storCounterX + 3;
    
end

for storEntries2 = 1:size(erp_story_pre,1)
    
    disp(strcat('Evaluating Story Post:', int2str(storEntries2)))
    
    my_data{storCounterA, 1} = erp_story_post{storEntries2,1}.erp_A_post;
    my_data{storCounterB, 1} = erp_story_post{storEntries2,1}.erp_B_post;
    my_data{storCounterX, 1} = erp_story_post{storEntries2,1}.erp_X_post;
    
    storCounterA = storCounterA + 3;
    storCounterB = storCounterB + 3;
    storCounterX = storCounterX + 3;
    
end


%% Time-lock analysis preparation

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

    qsubfeval('ERP_basecont_qsubRT', save_dir, subj, data_A_pre, data_A_post, data_B_pre, data_B_post, data_X_pre, data_X_post, 'memreq', 24*(1024^3), 'timreq',1*60*60,'memoverhead', 24*(1024^3), 'batchid',sprintf('TF_%s',subj)); 

    keep iSubject subjects root_dir save_dir
    
end

%% Load Data into structure

all_A_pre = {};
all_A_post = {};
for iSubject = 1:5%length(subjects)
    
    subj = subjects{iSubject};
    load (fullfile(save_dir, [subj,'_ERP_basecon']), 'erp_A_pre','erp_A_post');
    
    all_A_pre{iSubject}     = erp_A_pre;
    all_A_post{iSubject}    = erp_A_post;
    
    disp(int2str(iSubject))
    
    keep iSubject subjects root_dir save_dir all_A_pre all_A_post
    
    disp(strcat('Subj: ', int2str(iSubject)))
     
end

%% Searchlight stats

cfg             = [];
cfg.method      = 'triangulation';
cfg.layout      = 'CTF275_helmet.mat';
cfg.channel     = all_A_pre{1}.label;
neighbours      = ft_prepare_neighbours(cfg);

cfg = [] ;  
cfg.method      = 'mvpa';
cfg.searchlight = 'yes';
cfg.design      = [ones(length(all_A_pre{2}.trialinfo),1); 2*ones(length(all_A_post{2}.trialinfo),1)];
cfg.latency     = 'all';
cfg.avgovertime = 'yes';

cfg.mvpa.neighbours  = neighbours;

stat = ft_timelockstatistics(cfg, all_A_pre{2}, all_A_post{2})

%%

cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = 'CTF275_helmet.mat';            
cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
ft_topoplotER(cfg, stat);

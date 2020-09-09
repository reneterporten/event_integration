function [my_data] = rt_mytimelock_cov2(root_dir, subj)

load(fullfile(root_dir, [subj,'_dataclean.mat']));

% Recode trialinfo for linking events
indexLink   = 1:18:length(dataclean.trialinfo);
do_dont = -1;
for indices = 1:24  %Number of stories * 2
    if do_dont == 1
        dataclean.trialinfo(indexLink(indices):indexLink(indices)+5, 2) = 4;
        dataclean.trialinfo(indexLink(indices):indexLink(indices)+5, 3) = 4;
        dataclean.trialinfo(indexLink(indices):indexLink(indices)+5, 5) = 3;
        indexLink = indexLink+6;
    end
    do_dont = do_dont*-1;
end

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
    cfg         = [];
    cfg.trials  = find(data_post.trialinfo(:,4) == storyNr);
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


% Time-lock analysis of the phase x story x video x trial data
erp_story_pre   = {};
erp_story_post  = {};
for prepost = 1:size(data_pre_story_A,1)
    
    disp(strcat('Evaluating Story:', int2str(prepost)))
    
    
    cfg             = [];
    cfg.channel     = 'MEG';
    cfg.trials      = 'all';
    cfg.keeptrials  = 'yes';
    cfg.covariance  = 'yes';
    %cfg.preproc.demean          = 'yes';
    %cfg.preproc.baselinewindow  = [-0.3 0];
    erp_A_pre       = ft_timelockanalysis(cfg, data_pre_story_A{prepost});
    erp_B_pre       = ft_timelockanalysis(cfg, data_pre_story_B{prepost});
    erp_X_pre       = ft_timelockanalysis(cfg, data_pre_story_X{prepost});
    
    erp_A_post      = ft_timelockanalysis(cfg, data_pre_story_A{prepost});
    erp_B_post      = ft_timelockanalysis(cfg, data_pre_story_B{prepost});
    erp_X_post      = ft_timelockanalysis(cfg, data_pre_story_X{prepost});

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


% Arrange data structure for RSA
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

disp('Saving larger data structure...')
save(fullfile('/project/3012026.13/processed_RT/source reconstruction/', subj, 'tlck_structure.mat'), 'my_data', '-v7.3')
disp('Done.')

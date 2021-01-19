function [my_data] = rt_mytimelockv2(root_dir, subj)


% load data, remove bad trials and demean & detrend
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

cfg              = [];
cfg.feedback     = 'no';
cfg.method       = 'template';
cfg.template     = 'CTF275_neighb.mat';
cfg.planarmethod = 'sincos';
cfg.channel      = {'MEG'};
cfg.trials       = 'all';
cfg.neighbours   = ft_prepare_neighbours(cfg, data_pre);

% Disable planar combination
data_pre_tl      = ft_megplanar(cfg, data_pre);
data_post_tl     = ft_megplanar(cfg, data_post);

cfg              = [];
cfg.channel      = 'MEG';
cfg.trials       = 'all';
cfg.keeptrials   = 'yes';
cfg.lpfilter     = 'yes'; % apply lowpass filter
cfg.lpfreq       = 35;
data_pre_tl      = ft_timelockanalysis(cfg, data_pre_tl);
data_post_tl     = ft_timelockanalysis(cfg, data_post_tl);
% 
% cfg              = [];
% cfg.method       = 'sum';
% data_pre_tl      = ft_combineplanar(cfg,data_pre_tl);
% data_post_tl     = ft_combineplanar(cfg,data_post_tl);

cfg              = [];
cfg.latency      = [0 2.0];
data_pre_tl      = ft_selectdata(cfg,data_pre_tl);
data_post_tl     = ft_selectdata(cfg,data_post_tl);

% Run through all the stories
data_pre_story  = {};
data_post_story = {};
for storyNr = 1:max(unique(data_post.trialinfo(:,4)))
    
    cfg         = [];
    cfg.trials  = find(data_pre_tl.trialinfo(:,4) == storyNr);
    data_pre_story{storyNr,1} = ft_selectdata(cfg, data_pre_tl);
    
    % It appears that the linking event has been coded in as well (first 6 trials of post)
    cfg         = [];
    cfg.trials  = find(data_post_tl.trialinfo(:,4) == storyNr);
    data_post_story{storyNr,1} = ft_selectdata(cfg, data_post_tl);
    
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

% Arrange data structure for RSA
my_data     = {};
storCounterA = 1;
storCounterB = 2;
storCounterX = 3;
for storEntries = 1:size(data_pre_story_A,1)
    
    disp(strcat('Evaluating Story Pre:', int2str(storEntries)))
    
    my_data{storCounterA, 1} = data_pre_story_A{storEntries,1};
    my_data{storCounterB, 1} = data_pre_story_B{storEntries,1};
    my_data{storCounterX, 1} = data_pre_story_X{storEntries,1};
    
    storCounterA = storCounterA + 3;
    storCounterB = storCounterB + 3;
    storCounterX = storCounterX + 3;
    
end

for storEntries2 = 1:size(data_post_story_A,1)
    
    disp(strcat('Evaluating Story Post:', int2str(storEntries2)))
    
    my_data{storCounterA, 1} = data_post_story_A{storEntries2,1};
    my_data{storCounterB, 1} = data_post_story_B{storEntries2,1};
    my_data{storCounterX, 1} = data_post_story_X{storEntries2,1};
    
    storCounterA = storCounterA + 3;
    storCounterB = storCounterB + 3;
    storCounterX = storCounterX + 3;
    
end

clear data_pre_story_A data_pre_story_B data_pre_story_X
clear data_post_story_A data_post_story_B data_post_story_X 
clear dataclean data_pre data_post data_pre_tl data_post_tl
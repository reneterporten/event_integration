function [my_data] = rt_mytimelockv4(root_dir, subj, varargin)


% load data, remove bad trials and filter
load(fullfile(root_dir, [subj,'_dataclean.mat']), 'dataclean');

%%
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

%%
% do some preprocessing
cfg         = ft_getopt(varargin, 'cfg_preproc', []);
if ~any(contains(fieldnames(struct(cfg)), 'filter'))
  cfg.lpfilter      = 'yes';
  cfg.lpfreq        = 35;
  cfg.lpfilttype    = 'firws';
end
cfg.trials  = find(dataclean.trialinfo(:,5)==1);
data_pre    = ft_preprocessing(cfg, dataclean);

cfg.trials  = find(dataclean.trialinfo(:,5)==2);
data_post   = ft_preprocessing(cfg, dataclean);
clear dataclean;

cfg              = [];
cfg.feedback     = 'no';
cfg.method       = 'template';
cfg.template     = 'CTF275_neighb.mat';
cfg.planarmethod = 'sincos';
cfg.channel      = {'MEG'};
cfg.trials       = 'all';
cfg.neighbours   = ft_prepare_neighbours(cfg, data_pre);

% compute synthetic planar gradient
data_pre_tl      = ft_megplanar(cfg, data_pre);
data_post_tl     = ft_megplanar(cfg, data_post);
nstory           = max(unique(data_post.trialinfo(:,4)));
clear data_pre data_post;

% convert to a timelock representation
cfg              = ft_getopt(varargin, 'cfg_timelock', []);
cfg.channel      = 'MEG';
cfg.trials       = 'all';
cfg.keeptrials   = 'yes';
data_pre_tl      = ft_timelockanalysis(cfg, data_pre_tl);
data_post_tl     = ft_timelockanalysis(cfg, data_post_tl);

% prune the epochs
cfg              = [];
cfg.latency      = [-0.2 2]; %[0 2.0];
data_pre_tl      = ft_selectdata(cfg,data_pre_tl);
data_post_tl     = ft_selectdata(cfg,data_post_tl);

%%
% Run through all the stories, reorganise as a cell-array, one cell per
% story
data_pre_story  = cell(nstory, 1);
data_post_story = cell(nstory, 1);
for storyNr = 1:nstory
  
  cfg         = [];
  cfg.trials  = find(data_pre_tl.trialinfo(:,4) == storyNr);
  data_pre_story{storyNr,1} = ft_selectdata(cfg, data_pre_tl);
  
  % It appears that the linking event has been coded in as well (first 6 trials of post)
  cfg         = [];
  cfg.trials  = find(data_post_tl.trialinfo(:,4) == storyNr);
  data_post_story{storyNr,1} = ft_selectdata(cfg, data_post_tl);
  
end
clear data_pre_tl data_post_tl;

%%
% Reorganize the data into the individual events per story
data_pre_story_A = cell(nstory,1);
data_pre_story_B = cell(nstory,1);
data_pre_story_X = cell(nstory,1);

data_post_story_A = cell(nstory,1);
data_post_story_B = cell(nstory,1);
data_post_story_X = cell(nstory,1);

for storyNr = 1:nstory
  
  % Pre data
  cfg                         = [];
  cfg.trials                  = find(data_pre_story{storyNr,1}.trialinfo(:,3) == 1);
  data_pre_story_A{storyNr,1} = ft_selectdata(cfg, data_pre_story{storyNr,1});
  
  cfg                         = [];
  cfg.trials                  = find(data_pre_story{storyNr,1}.trialinfo(:,3) == 2);
  data_pre_story_B{storyNr,1} = ft_selectdata(cfg, data_pre_story{storyNr,1});
  
  cfg                         = [];
  cfg.trials                  = find(data_pre_story{storyNr,1}.trialinfo(:,3) == 3);
  data_pre_story_X{storyNr,1} = ft_selectdata(cfg, data_pre_story{storyNr,1});
  
  % Post data
  cfg                             = [];
  cfg.trials                      = find(data_post_story{storyNr,1}.trialinfo(:,3) == 1);
  data_post_story_A{storyNr,1}    = ft_selectdata(cfg, data_post_story{storyNr,1});
  
  cfg                             = [];
  cfg.trials                      = find(data_post_story{storyNr,1}.trialinfo(:,3) == 2);
  data_post_story_B{storyNr,1}    = ft_selectdata(cfg, data_post_story{storyNr,1});
  
  cfg                             = [];
  cfg.trials                      = find(data_post_story{storyNr,1}.trialinfo(:,3) == 3);
  data_post_story_X{storyNr,1}    = ft_selectdata(cfg, data_post_story{storyNr,1});
  
end
clear data_pre_story data_post_story;

%%
% Arrange data structure for RSA
my_data      = cell(nstory*6, 1);
storCounterA = 1;
storCounterB = 2;
storCounterX = 3;
for storEntries = 1:nstory
  
  disp(strcat('Evaluating Story Pre:', int2str(storEntries)))
  
  my_data{storCounterA, 1} = data_pre_story_A{storEntries,1};
  my_data{storCounterB, 1} = data_pre_story_B{storEntries,1};
  my_data{storCounterX, 1} = data_pre_story_X{storEntries,1};
  
  storCounterA = storCounterA + 3;
  storCounterB = storCounterB + 3;
  storCounterX = storCounterX + 3;
  
end

for storEntries2 = 1:nstory
  
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

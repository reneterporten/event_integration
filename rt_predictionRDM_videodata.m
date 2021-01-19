function [videoDataRSM_magn_conc, videoDataRSM_phas_conc] = rt_predictionRDM_videodata(data)

% Load Data
root_dir    = '/project/3012026.13/processed/';
script_dir  = '/project/3012026.13/scripts_RT/';
dataclean   = data;
load(fullfile(script_dir, 'video_data.mat'));

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

% Select trials for pre and post
cfg         = [];
cfg.trials  = find(dataclean.trialinfo(:,5)==1);
data_pre    = ft_selectdata(cfg, dataclean);

cfg         = [];
cfg.trials  = find(dataclean.trialinfo(:,5)==2);
data_post   = ft_selectdata(cfg, dataclean);


% % Check whehter trialstructure has a length of 216, if not the linking
% % events have been added to the structure
% if length(data_pre.trialinfo) > 216
%     
%     disp('Adjusting trialfino in pre-data, more than 216 trials detected.')
%     selectpre = true(length(data_pre.trialinfo),1);
%     adjustpre = [false(6,1); true(18,1)];
%     
%     for predata = 1:24:length(data_pre.trialinfo)
%         selectpre(predata:predata+23,1) = adjustpre;
%     end
%     
%     cfg         = [];
%     cfg.trials  = selectpre;
%     data_pre   = ft_selectdata(cfg, data_pre);
%     
% end
% 
% if length(data_post.trialinfo) > 216
%     
%     disp('Adjusting trialfino in post-data, more than 216 trials detected.')
%     selectpost = true(length(data_post.trialinfo),1);
%     adjustpost = [false(6,1); true(18,1)];
%     
%     for postdata = 1:24:length(selectpost)
%         selectpost(postdata:postdata+23,1) = adjustpost;
%     end
%     
%     cfg         = [];
%     cfg.trials  = selectpost;
%     data_post   = ft_selectdata(cfg, data_post);
%     
% end
% 
% 
% Get video data for pre trials
trials_magn_pre = [];
trials_phas_pre = [];
for trials = 1:size(data_pre.trialinfo, 1)
    
    % Loop through trials and create video indices to be read from the
    % video data
    storyinfo   = data_pre.trialinfo(trials, 4);
    videoinfo   = data_pre.trialinfo(trials, 3);
    
    storyidx        = video_data.storyidx == storyinfo;
    magn_data_story = video_data.magnitude(storyidx,:);
    phas_data_story = video_data.phase(storyidx,:);
    trials_magn_pre = [trials_magn_pre; magn_data_story(videoinfo,:)];
    trials_phas_pre = [trials_phas_pre; phas_data_story(videoinfo,:)];
    
end

% Get vide data for post trials
trials_magn_post = [];
trials_phas_post = [];
for trials = 1:size(data_pre.trialinfo, 1)
    
    % Loop through trials and create video indices to be read from the
    % video data
    storyinfo   = data_post.trialinfo(trials, 4);
    videoinfo   = data_post.trialinfo(trials, 3);
    
    storyidx        = video_data.storyidx == storyinfo;
    magn_data_story = video_data.magnitude(storyidx,:);
    phas_data_story = video_data.phase(storyidx,:);
    trials_magn_post = [trials_magn_post; magn_data_story(videoinfo,:)];
    trials_phas_post = [trials_phas_post; phas_data_story(videoinfo,:)];
    
end

% Put trialinfo and video data into one data structure for ordering
% Trial, video event, story, all phase/magnitude values per frame
pre_trialinfo_magn  = [data_pre.trialinfo(:,1), data_pre.trialinfo(:,3:4), trials_magn_pre];
pre_trialinfo_phas  = [data_pre.trialinfo(:,1), data_pre.trialinfo(:,3:4), trials_phas_pre];
post_trialinfo_magn = [data_post.trialinfo(:,1), data_post.trialinfo(:,3:4), trials_magn_post];
post_trialinfo_phas = [data_post.trialinfo(:,1), data_post.trialinfo(:,3:4), trials_phas_post];

% Sort the trialinfo structure by story and then by video event
pre_trialinfo_magn  = sortrows(pre_trialinfo_magn, [3, 2, 1]);
pre_trialinfo_phas  = sortrows(pre_trialinfo_phas, [3, 2, 1]);
post_trialinfo_magn = sortrows(post_trialinfo_magn, [3, 2, 1]);
post_trialinfo_phas = sortrows(post_trialinfo_phas, [3, 2, 1]);

% Get rid of the additional trial info and only leave video data
pre_trialinfo_magn  = pre_trialinfo_magn(:,4:end);
pre_trialinfo_phas  = pre_trialinfo_phas(:,4:end);
post_trialinfo_magn = post_trialinfo_magn(:,4:end);
post_trialinfo_phas = post_trialinfo_phas(:,4:end);

% Calculate the variance over frames
pre_trialinfo_magn_col  = var(pre_trialinfo_magn');
pre_trialinfo_phas_col  = var(pre_trialinfo_phas');
post_trialinfo_magn_col = var(post_trialinfo_magn');
post_trialinfo_phas_col = var(post_trialinfo_phas');

pre_trialinfo_magn_row  = pre_trialinfo_magn_col';
pre_trialinfo_phas_row  = pre_trialinfo_phas_col';
post_trialinfo_magn_row = post_trialinfo_magn_col';
post_trialinfo_phas_row = post_trialinfo_phas_col';

% Concatonate magn and phase matrices for pre and post
trialinfo_magn_col = [pre_trialinfo_magn_col, post_trialinfo_magn_col];
trialinfo_phas_col = [pre_trialinfo_phas_col, post_trialinfo_phas_col];
trialinfo_magn_row = [pre_trialinfo_magn_row; post_trialinfo_magn_row];
trialinfo_phas_row = [pre_trialinfo_phas_row; post_trialinfo_phas_row];

% Create big prediction matrix
videoDataRSM_magn_conc = nan(length(trialinfo_magn_row), length(trialinfo_magn_col));
videoDataRSM_phas_conc = nan(length(trialinfo_magn_row), length(trialinfo_magn_col));

% Include across phase comparisons
% Loop through pre matrix for magn and phase
for matRow = 1:1:length(videoDataRSM_magn_conc)
    for matCol = 1:1:length(videoDataRSM_magn_conc)
        
        videoDataRSM_magn_conc(matRow, matCol) = 1-(trialinfo_magn_row(matRow)-trialinfo_magn_col(matCol));
        videoDataRSM_phas_conc(matRow, matCol) = 1-(trialinfo_phas_row(matRow)-trialinfo_phas_col(matCol));
        
    end
end
% 
% % Do NOT include across phase comparisons
% % Big matrix pre
% big_magn_pre = nan(18*12, 18*12);
% big_phas_pre = nan(18*12, 18*12);
% 
% % Big matrix post
% big_magn_post = nan(18*12, 18*12);
% big_phas_post = nan(18*12, 18*12);
% 
% % Loop through pre matrix for magn and phase
% for preMatRow = 1:1:length(big_matrix_pre)
%     for preMatCol = 1:1:length(big_matrix_pre)
%         
%         big_magn_pre(preMatRow, preMatCol) = 1-(pre_trialinfo_magn_row(preMatRow)-pre_trialinfo_magn_col(preMatCol));
%         big_phas_pre(preMatRow, preMatCol) = 1-(pre_trialinfo_phas_row(preMatRow)-pre_trialinfo_phas_col(preMatCol));
%         
%     end
% end
% 
% % Loop through post matrix
% for postMatRow = 1:1:length(big_matrix_post)
%     for postMatCol = 1:1:length(big_matrix_post)
%         
%         big_magn_post(postMatRow, postMatCol) = 1-(post_trialinfo_magn_row(postMatRow)-post_trialinfo_magn_col(postMatCol));
%         big_phas_post(postMatRow, postMatCol) = 1-(post_trialinfo_phas_row(postMatRow)-post_trialinfo_phas_col(postMatCol));
%         
%     end
% end
% 
% % Create big prediction matrix
% lengthData = length(big_magn_pre);
% videoDataRSM_magn = nan(length(big_magn_pre)*2, length(big_magn_pre)*2);
% videoDataRSM_phas = nan(length(big_magn_pre)*2, length(big_magn_pre)*2);
% 
% videoDataRSM_magn(1:lengthData, 1:lengthData)           = big_magn_pre;
% videoDataRSM_magn(lengthData+1:end, lengthData+1:end)   = big_magn_post;

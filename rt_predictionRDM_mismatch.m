function [mismatchRSM] = rt_predictionRDM_mismatch(data)

% Load Data
root_dir    = '/project/3012026.13/processed/';
dataclean   = data;

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

% Loop through trials and create time-distance structure
pre_distance            = nan(size(data_pre.trialinfo, 1),1);
post_distance           = nan(size(data_post.trialinfo, 1),1);
pre_distance_counter    = 0;
pre_counter_start       = 0;
post_distance_counter   = 0;
post_counter_start      = 0;
pre_story               = data_pre.trialinfo(1,4);
post_story              = data_post.trialinfo(1,4);
for trials = 1:size(data_pre.trialinfo, 1)
    
    % Check whether a new story began, if so set start counter back to 0
    if data_pre.trialinfo(trials,4) ~= pre_story
        pre_counter_start       = 0;
        pre_distance_counter    = 0;
        pre_story               = data_pre.trialinfo(trials,4);
    end
    
    if data_post.trialinfo(trials,4) ~= post_story
        post_counter_start      = 0;
        post_distance_counter   = 0;
        post_story              = data_post.trialinfo(trials,4);
    end
    
    % Check when first A trial occurs in pre
    if data_pre.trialinfo(trials,3) == 1
        pre_counter_start = 1;
    else
        pre_distance(trials,1) = 4; % Max distance for first trials
    end
    
    % Check when first A trial occurs in post
    if data_post.trialinfo(trials,3) == 1
        post_counter_start = 1;
    else
        post_distance(trials,1) = 4; % Max distance for first trials
    end
    
    % Apply distance scores in pre
    if pre_counter_start == 1 && data_pre.trialinfo(trials,3) == 1
        pre_distance(trials,1) = 0;
        pre_distance_counter   = 0;
    elseif pre_counter_start == 1
        pre_distance_counter   = pre_distance_counter + 2;
        pre_distance(trials,1) = pre_distance_counter;
    end
    
    % Apply distance scores in post
    if post_counter_start == 1 && data_post.trialinfo(trials,3) == 1
        post_distance(trials,1) = 0;
        post_distance_counter   = 0;
    elseif post_counter_start == 1
        post_distance_counter   = post_distance_counter + 2;
        post_distance(trials,1) = post_distance_counter;
    end
    
end

% Applay distance structure to trialinfo structure of pre and post
data_pre.trialinfo  = [data_pre.trialinfo, pre_distance];
data_post.trialinfo = [data_post.trialinfo, post_distance];

% Loop through the data per story and sort by A, B & X
sorted_pre  = [];
sorted_post = [];
for currentStory = 1:12
    
    idxStoryPre     = find(data_pre.trialinfo(:,4) == currentStory);
    idxStoryPost    = find(data_post.trialinfo(:,4) == currentStory);
    sorted_pre      = [sorted_pre; sortrows(data_pre.trialinfo(idxStoryPre,:), 3)];
    sorted_post     = [sorted_post; sortrows(data_post.trialinfo(idxStoryPost,:), 3)];
    
end

% Raw Nan
nan_matrix = nan(18,18);

% Pre Matrix
pre_matrix              = nan_matrix;
pre_matrix(1:6, 7:12)   = -1;
pre_matrix(1:6, 13:18)  = 1;
pre_matrix(7:12, 1:6)   = -1;
pre_matrix(13:18, 1:6)  = 1;

% Post Matrix
post_matrix = pre_matrix.*-1;

% Big matrix pre
big_matrix_pre = nan(18*12, 18*12);

% Big matrix post
big_matrix_post = nan(18*12, 18*12);

% Loop through pre matrix
for preMat = 1:18:length(big_matrix_pre)

    big_matrix_pre(preMat:preMat+17, preMat:preMat+17) = pre_matrix ;
    
end

% Loop through post matrix
for postMat = 1:18:length(big_matrix_post)
    
    big_matrix_post(postMat:postMat+17, postMat:postMat+17) = post_matrix ;
    
end

% With the raw matrices for pre and post being ready, absolute differences 
% in time can be calculated and stored in the specific cells
% Pre Phase
sorted_pre_row = sorted_pre(:,7);
sorted_pre_col = sorted_pre_row';
for prerows = 1:length(big_matrix_pre)
    for precols = 1:length(big_matrix_pre)
        
        big_matrix_pre(prerows, precols) = big_matrix_pre(prerows, precols)*abs(sorted_pre_row(prerows)-sorted_pre_col(precols));
        
    end
end

% Post Phase
sorted_post_row = sorted_post(:,7);
sorted_post_col = sorted_post_row';
for postrows = 1:length(big_matrix_post)
    for postcols = 1:length(big_matrix_post)
        
        big_matrix_post(postrows, postcols) = big_matrix_post(postrows, postcols)*abs(sorted_post_row(postrows)-sorted_post_col(postcols));
        
    end
end

% Combine pre and post phase matrices into one big prediction RSM
leng_data = length(big_matrix_pre);
mismatchRSM                                                     = nan(length(big_matrix_pre)*2, length(big_matrix_post)*2);
mismatchRSM(1:leng_data, 1:leng_data)                           = big_matrix_pre;
mismatchRSM(leng_data+1:leng_data*2, leng_data+1:leng_data*2)   = big_matrix_post;

%imagesc(mismatchRSM)

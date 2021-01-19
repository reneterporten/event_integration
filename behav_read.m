%% Read in behavioral rating data
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

%%
behav_data = readtable('/project/3012026.13/processed_RT/logfiles/sub004_Relatedness_answers.txt');
behav_data = table2array(behav_data(:,1:6));

% Sort data by story
[B, I]          = sort(behav_data,1);
sortStory_data  = behav_data(I(:,1),:);

% Identify pre and post data
pre_data_idx    = sortStory_data(:,2) == 1;
post_data_idx   = sortStory_data(:,2) == 3;

pre_data        = sortStory_data(pre_data_idx,:);
post_data       = sortStory_data(post_data_idx,:);

for preStory = 1:12
    % Prepare matrix for comparison
    pre_mat = ones(3,3);  

    pre_pairs = int2str(pre_data(:,3));

    pre_story_idx   = pre_data(:,1) == preStory;
    pre_story_data  = pre_data(pre_story_idx,:);

    for preStorNum = 1:size(pre_story_data,1)

        idx_4_mat_str   = int2str(pre_story_data(preStorNum, 3));
        idx_4_mat       = [str2num(idx_4_mat_str(1,1)), str2num(idx_4_mat_str(1,2))];
        preStorRating   = pre_story_data(preStorNum,6);

        pre_mat(idx_4_mat(1,1), idx_4_mat(1,2)) = preStorRating;

    end

    % Calculate mean rating between comparisons
    low_pre_mat     = tril(pre_mat, -1);
    upp_pre_mat     = triu(pre_mat, 1)';
    mean_pre_mat    = (low_pre_mat + upp_pre_mat)./2;
    mean_pre_mat    = mean_pre_mat + mean_pre_mat';
    diagonalidx     = logical(eye(size(mean_pre_mat)));
    mean_pre_mat(diagonalidx) = 6;

    % Create big predictions matrix for pre and post
    trial_pre_mat = nan(18,18);
    trial_mad_idx = 1:6:18;
    for rowspre = 1:size(mean_pre_mat,1)

        all_A = ones(6,6)*mean_pre_mat(rowspre,1);
        all_B = ones(6,6)*mean_pre_mat(rowspre,2);
        all_X = ones(6,6)*mean_pre_mat(rowspre,3);
        all_cond = [all_A, all_B, all_X];

        trial_pre_mat(trial_mad_idx(rowspre):trial_mad_idx(rowspre)+5,:) = all_cond;

    end
    
    all_data_pre{preStory} = trial_pre_mat;

end

for postStory = 1:12
    % Prepare matrix for comparison
    post_mat = ones(3,3);  

    post_pairs = int2str(post_data(:,3));

    post_story_idx   = post_data(:,1) == postStory;
    post_story_data  = post_data(post_story_idx,:);

    for postStorNum = 1:size(post_story_data,1)

        idx_4_mat_str   = int2str(post_story_data(postStorNum, 3));
        idx_4_mat       = [str2num(idx_4_mat_str(1,1)), str2num(idx_4_mat_str(1,2))];
        postStorRating   = post_story_data(postStorNum,6);

        post_mat(idx_4_mat(1,1), idx_4_mat(1,2)) = postStorRating;

    end

    % Calculate mean rating between comparisons
    low_post_mat     = tril(post_mat, -1);
    upp_post_mat     = triu(post_mat, 1)';
    mean_post_mat    = (low_post_mat + upp_post_mat)./2;
    mean_post_mat    = mean_post_mat + mean_post_mat';
    diagonalidx      = logical(eye(size(mean_post_mat)));
    mean_post_mat(diagonalidx) = 6;

    % Create big predictions matrix for pre and post
    trial_post_mat = nan(18,18);
    trial_mad_idx = 1:6:18;
    for rowspost = 1:size(mean_post_mat,1)

        all_A = ones(6,6)*mean_post_mat(rowspost,1);
        all_B = ones(6,6)*mean_post_mat(rowspost,2);
        all_X = ones(6,6)*mean_post_mat(rowspost,3);
        all_cond = [all_A, all_B, all_X];

        trial_post_mat(trial_mad_idx(rowspost):trial_mad_idx(rowspost)+5,:) = all_cond;

    end
    
    all_data_post{postStory} = trial_post_mat;

end


all_behav_data  = nan(12*18*2,12*18*2);
story_counter   = 1;
for preRDMs = 1:18:(size(all_behav_data,1)/2)
    
    all_behav_data(preRDMs:(preRDMs+17), preRDMs:(preRDMs+17)) = all_data_pre{story_counter};
    story_counter = story_counter + 1;
    
end

story_counter   = 1;
for postRDMs = 217:18:(size(all_behav_data,1))
    
    all_behav_data(postRDMs:(postRDMs+17), postRDMs:(postRDMs+17)) = all_data_post{story_counter};
    story_counter = story_counter + 1;
    
end

imagesc(all_behav_data)

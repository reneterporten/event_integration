function [link_mech] = rt_predictionRDM_linkmech()

% Raw Nan
nan_matrix = nan(18,18);

% Pre Matrix
pre_matrix              = nan_matrix;
pre_matrix(1:6, 7:12)   = -1;
pre_matrix(1:6, 13:18)  = 1;
pre_matrix(7:12, 1:6)   = -1;
pre_matrix(7:12, 13:18) = 1;
pre_matrix(13:18, 1:6)  = 1;
pre_matrix(13:18, 7:12) = 1;

% Post Matrix
post_matrix = pre_matrix.*-1;

% Big matrix pre
big_matrix_pre = nan(18*12, 18*12);

% Big matrix post
big_matrix_post = nan(18*12, 18*12);

% Loop through pre matrix
for preMatRow = 1:18:length(big_matrix_pre)
    
    for preMatCol = 1:18:length(big_matrix_pre)
        big_matrix_pre(preMatRow:preMatRow+17, preMatCol:preMatCol+17) = pre_matrix ;
    end
    
end

% Loop through post matrix
for postMatRow = 1:18:length(big_matrix_post)
    
    for postMatCol = 1:18:length(big_matrix_post)
        big_matrix_post(postMatRow:postMatRow+17, postMatCol:postMatCol+17) = post_matrix ;
    end
    
end

% Loop through pre and post and insert within story NaN
for nanTrain =  1:18:length(big_matrix_pre)
   
    big_matrix_pre(nanTrain:nanTrain+17, nanTrain:nanTrain+17) = nan_matrix;
    big_matrix_post(nanTrain:nanTrain+17, nanTrain:nanTrain+17) = nan_matrix;
    
end

% Prediction matrix
leng_data = length(big_matrix_pre);

link_mech = nan(length(big_matrix_pre)*2, length(big_matrix_post)*2);
link_mech(1:leng_data, 1:leng_data) = big_matrix_pre;
link_mech(leng_data+1:leng_data*2, leng_data+1:leng_data*2) = big_matrix_post;


%link_mech(isnan(link_mech)) = 0;
%imagesc(link_mech)






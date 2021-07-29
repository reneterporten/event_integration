function rt_corrbehav(varargin)

% Function to correlate behavioral relatedness ratings with imaging data

saveflag    = ft_getopt(varargin, 'saveflag', false);
savepath    = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
datadir     = ft_getopt(varargin, 'datadir', '/project/3012026.13/logfiles/MEG/');
savename    = ft_getopt(varargin, 'savename', 'coherence'); 
suff        = ft_getopt(varargin, 'suff', 'Relatedness_answers.txt1');
suffcoh     = ft_getopt(varargin, 'suffcoh', 'coh.mat');
atlasrois   = ft_getopt(varargin, 'atlasrois', {'A10m', 'A11m', 'A13', 'A14m', 'A32sg', 'Hipp'}); % Either 'all' or cell array with ROIs
fname       = ft_getopt(varargin, 'fname', 'cohspctrm');

cd(datadir);

% Organize behavioral data
pairs = [12, 21; 13, 31; 23, 32]; % AB; AX; BX
d = dir(sprintf('sub*%s*',suff));
for k = 1:numel(d)  
    behav_data = readtable(d(k).name, 'FileType', 'text');
    behav_data = removevars(behav_data,{'Var7'});
    behav_data = table2array(behav_data);
    pre_idx     = behav_data(:,2) == 1;
    post_idx    = behav_data(:,2) == 3;
    behav_pre   = behav_data(pre_idx,:);
    behav_post  = behav_data(post_idx,:);
    for paircount = 1:size(pairs, 1)
        pre_comp(paircount,k)  = mean(behav_pre(ismember(behav_pre(:,3), pairs(paircount, :)), 6));
        post_comp(paircount,k) = mean(behav_post(ismember(behav_post(:,3), pairs(paircount, :)), 6));
    end
end

% Loading and organizing imaging data
cd(savepath)
d = dir(sprintf('sub*%s*',suffcoh));
for k = 1:numel(d)
    connectivity_data = load(d(k).name);
    sorted_idx = get_sortedrois(connectivity_data.coh{1,1}.label, atlasrois);    
    for i = 1:length(connectivity_data.conlabel)
        current_coh = connectivity_data.coh{1,i};
        current_coh.(fname) = current_coh.(fname)(sorted_idx, sorted_idx);
        P   = blkdiag(ones(1,5)./5,ones(1,2)./2,ones(1,2)./2,ones(1,5)./5);
        tmp = current_coh.(fname);
        tmp = P*tmp*P';
        tmp = tmp - diag(diag(tmp));
        coh_data{k,i}.(fname) = squareform(tmp).';
        coh_data{k,i}.dimord  = 'chan_freq';
        coh_data{k,i}.grad    = connectivity_data.coh{1,i}.grad; 
        coh_data{k,i}.freq    = connectivity_data.coh{1,i}.freq;
        coh_data{k,i}.label   = {'Hippl_MPFCl';'Hippr_MPFCl';'MPFCr_MPFCl';'Hippr_Hippl';'MPFCr_Hippl';'Hippr_MPFCr'};
    end
end

% Average coherence data to create pairs comparable to behavioral data
for subs = 1:size(coh_data,1)
    for cond = 1:size(coh_data,2)
        coh_pairs{subs, cond} = coh_data{subs,cond};
    end
    % Change this into collapsing conditions into one 6x1 fields as is the
    % behavioral data
    % Pre
    coh_pairs{subs, 1}.cohspctrm = (coh_data{subs,1}.cohspctrm + coh_data{subs,2}.cohspctrm)/2; % AB
    coh_pairs{subs, 2}.cohspctrm = (coh_data{subs,1}.cohspctrm + coh_data{subs,3}.cohspctrm)/2; % AX
    coh_pairs{subs, 3}.cohspctrm = (coh_data{subs,2}.cohspctrm + coh_data{subs,3}.cohspctrm)/2; % BX
    % Post
    coh_pairs{subs, 4}.cohspctrm = (coh_data{subs,4}.cohspctrm + coh_data{subs,5}.cohspctrm)/2; % AB
    coh_pairs{subs, 5}.cohspctrm = (coh_data{subs,4}.cohspctrm + coh_data{subs,6}.cohspctrm)/2; % AX
    coh_pairs{subs, 6}.cohspctrm = (coh_data{subs,5}.cohspctrm + coh_data{subs,6}.cohspctrm)/2; % BX
    for cond = 1:size(coh_data,2)
        coh_pairs{subs, cond}.cohspctrm_left = coh_pairs{subs, cond}.cohspctrm(1,1); % Hippl_MPFCl
        coh_pairs{subs, cond}.cohspctrm_right = coh_pairs{subs, cond}.cohspctrm(6,1); % Hippr_MPFCr
        coh_pairs{subs, cond}.relatedness = [pre_comp(:, subs); post_comp(:, subs)];
    end
end

% Correlate relatedness with coherence across subjects

% To-Do




disp('hi')



%% -------------------------- SUB FUNCTIONS -------------------------- %%

function [sorted_idx] = get_sortedrois(datalabels, atlasrois)

% Sorts the ROIs or all parcels to left and right hemisphere
% Provides index for further analyses

idx         = [];
leftidx     = [];
rightidx    = [];
if strcmp(atlasrois, 'all')
    atlasrois = datalabels;
end
for lab = 1:numel(datalabels)
    for roirun = 1:numel(atlasrois)
        if contains(datalabels{lab}, atlasrois(roirun))
            idx = [idx; lab];
            if contains(datalabels{lab}, 'Left')
                leftidx = [leftidx; true];
                rightidx = [rightidx; false];
            elseif contains(datalabels{lab}, 'Right')
                leftidx = [leftidx; false];
                rightidx = [rightidx; true];
            end
        end
    end
end
left_side   = idx(logical(leftidx));
right_side  = idx(logical(rightidx));
sorted_idx  = [left_side; flip(right_side)];




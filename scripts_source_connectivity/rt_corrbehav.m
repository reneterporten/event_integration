function rt_corrbehav(varargin)

% Function to correlate behavioral relatedness ratings with imaging data

saveflag    = ft_getopt(varargin, 'saveflag', false);
savepath    = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
datadir     = ft_getopt(varargin, 'datadir', '/project/3012026.13/logfiles/Relatedness/');
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
        pre_comp(k,paircount)  = mean(behav_pre(ismember(behav_pre(:,3), pairs(paircount, :)), 6));
        post_comp(k,paircount) = mean(behav_post(ismember(behav_post(:,3), pairs(paircount, :)), 6));
    end
    % Resort relatedness data for odd subjects as B-X switches
    if mod(str2double(extractBetween(d(k).name,'sub','_Relatedness')), 2) % in case of odd subject
        pre_comp(k,:) = pre_comp(k,[2, 1, 3]);
        post_comp(k,:) = post_comp(k,[2, 1, 3]);
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
    coh_pairs{subs, 1} = coh_data{subs, 1};
    coh_pairs{subs, 1}.cohpairs = zeros(6,6);
    % Pre
    coh_pairs{subs, 1}.cohpairs(:,1) = (coh_data{subs,1}.cohspctrm + coh_data{subs,2}.cohspctrm)/2; % AB
    coh_pairs{subs, 1}.cohpairs(:,2) = (coh_data{subs,1}.cohspctrm + coh_data{subs,3}.cohspctrm)/2; % AX
    coh_pairs{subs, 1}.cohpairs(:,3) = (coh_data{subs,2}.cohspctrm + coh_data{subs,3}.cohspctrm)/2; % BX
    % Post
    coh_pairs{subs, 1}.cohpairs(:,4) = (coh_data{subs,4}.cohspctrm + coh_data{subs,5}.cohspctrm)/2; % AB
    coh_pairs{subs, 1}.cohpairs(:,5) = (coh_data{subs,4}.cohspctrm + coh_data{subs,6}.cohspctrm)/2; % AX
    coh_pairs{subs, 1}.cohpairs(:,6) = (coh_data{subs,5}.cohspctrm + coh_data{subs,6}.cohspctrm)/2; % BX

    coh_pairs{subs, 1}.relatedness = [pre_comp(subs,:), post_comp(subs,:)];
    coh_pairs{subs, 1}.pair_dimord = 'chan_pair';
    coh_pairs{subs, 1}.pairlabel = {'AB_pre', 'AX_pre', 'BX_pre', 'AB_post', 'AX_post', 'BX_post' };
end

% Correlate relatedness with coherence across subjects
current_data_matrix = [];
for chan = 1:size(coh_pairs{1,1}.cohpairs, 1)
    for subs = 1:size(coh_pairs,1)
        current_data_matrix = [current_data_matrix; coh_pairs{subs,1}.cohpairs(chan,:), coh_pairs{subs,1}.relatedness];
    end
    data_matrix{chan,1} = current_data_matrix;
    current_data_matrix = [];
end

corr_matrix = zeros(6,6);
for chan = 1:6
    for pair = 1:6
        corr_matrix(chan,pair) = corr(data_matrix{chan}(:,pair), data_matrix{chan}(:,pair+6));
    end
end

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




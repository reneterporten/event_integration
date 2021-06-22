%% Pipeline to initiate RSA on source reconstructed data
%%

% Common fieltrip
%addpath /home/common/matlab/fieldtrip/
%addpath /home/common/matlab/fieldtrip/qsub/
%addpath /home/common/matlab/fieldtrip/private/
% JM forked fieldtrip
addpath /project/3012026.13/scripts_RT/git_event_integration/fieldtrip/
% Private scripts
addpath /project/3012026.13/scripts_RT/scripts_source_connectivity/
%addpath /project/3012026.13/scripts_RT/
ft_defaults

root_dir     = '/project/3012026.13/processed/';
save_dir     = '/project/3012026.13/processed_RT/source reconstruction/';

subjects = {'sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-022','sub-023',...
    'sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-031'...
    'sub-032','sub-033','sub-034','sub-035','sub-036','sub-037','sub-038'};


%% Volume realign

for subj = 1:length(subjects)
    
    disp(strcat('Subject aligning:', int2str(subj)))
    rt_volumealign(subjects{subj});
    disp(strcat('Subject aligning done:', int2str(subj)))
    
end


%% Create head- and sourcemodel

for subj = 1:length(subjects)
    
    disp(strcat('Subject headsource:', int2str(subj)))
    rt_headsource(subjects{subj});
    disp(strcat('Subject headsource done:', int2str(subj)))
    
end


%% Do source level anaylsis and estimate coherence

for subj = 1:length(subjects)
    
    disp(strcat('Subject sourceanalysis:', int2str(subj)))
    rt_sourcelevelanalysis(subjects{subj}, 'saveflag', true)
    disp(strcat('Subject sourceanalysis done:', int2str(subj)))
    
end


%% Estimate coherence

suff = '_coh.mat';
rt_collectsourcedata(suff, 'saveflag', true, 'connectivity', 'coh')
rt_collectsourcedata(suff, 'saveflag', true, 'connectivity', 'imcoh')
rt_collectsourcedata(suff, 'saveflag', true, 'connectivity', 'mim')


%% Plot coherence pre vs post, HC & mPFC

load(fullfile('/project/3012026.13/jansch/', 'groupdata_coh_coherence.mat'))
%load(fullfile('/project/3012026.13/jansch/', 'groupdata_imcoh_coherence.mat'))
%load(fullfile('/project/3012026.13/jansch/', 'groupdata_mim_coherence.mat'))

% Get index of ROIs, sorted by left and right hemisphere
idx = [];
leftidx = [];
rightidx = [];
for lab = 1:numel(Fcon{1}.label)
    if contains(Fcon{1}.label{lab}, 'MFG') || contains(Fcon{1}.label{lab}, 'Hipp')
        idx = [idx; lab];
        if contains(Fcon{1}.label{lab}, 'Left')
            leftidx = [leftidx; true];
            rightidx = [rightidx; false];
        elseif contains(Fcon{1}.label{lab}, 'Right')
            leftidx = [leftidx; false];
            rightidx = [rightidx; true];
        end
    end
end

left_side   = idx(logical(leftidx));
right_side  = idx(logical(rightidx));
sorted_idx  = [left_side; flip(right_side)];

%
Apre = Fcon{1}.cohspctrm(sorted_idx,:);
Apre = abs(Apre(:,sorted_idx));

Apost = Fcon{4}.cohspctrm(sorted_idx,:);
Apost = abs(Apost(:,sorted_idx));

imagesc(Apre, [0 0.4])
figure;
imagesc(Apost, [0 0.4])


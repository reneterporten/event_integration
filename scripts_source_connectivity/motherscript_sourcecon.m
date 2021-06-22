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

load(fullfile('/project/3012026.13/jansch/', 'groupdata_coh_coherence.mat')); fname = 'cohspctrm';
%load(fullfile('/project/3012026.13/jansch/', 'groupdata_imcoh_coherence.mat')); fname = 'cohspctrm';
%load(fullfile('/project/3012026.13/jansch/', 'groupdata_mim_coherence.mat')); fname = 'mimspctrm';

% Regions that fall under the MPFC
% A10m, A11m, A13, A14m, A32sg 
% including left (_l) and right (_r) variants

% Get index of ROIs, sorted by left and right hemisphere
idx = [];
leftidx = [];
rightidx = [];
roi_names = {'A10m', 'A11m', 'A13', 'A14m', 'A32sg', 'Hipp'};
for lab = 1:numel(Fcon{1}.label)
    for roirun = 1:numel(roi_names)
        if contains(Fcon{1}.label{lab}, roi_names(roirun))
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
end

left_side   = idx(logical(leftidx));
right_side  = idx(logical(rightidx));
sorted_idx  = [left_side; flip(right_side)];


%% Select condition specific data

Apre = Fcon{1}.(fname)(sorted_idx,:);
Apre = abs(Apre(:,sorted_idx));
Apost = Fcon{4}.(fname)(sorted_idx,:);
Apost = abs(Apost(:,sorted_idx));

Bpre = Fcon{2}.(fname)(sorted_idx,:);
Bpre = abs(Bpre(:,sorted_idx));
Bpost = Fcon{5}.(fname)(sorted_idx,:);
Bpost = abs(Bpost(:,sorted_idx));

Xpre = Fcon{3}.(fname)(sorted_idx,:);
Xpre = abs(Xpre(:,sorted_idx));
Xpost = Fcon{6}.(fname)(sorted_idx,:);
Xpost = abs(Xpost(:,sorted_idx));


%% Plot figures

figure;
% A
subplot(3,2,1);
imagesc(Apre); cax1=caxis;%, [0 0.4])
subplot(3,2,2);
imagesc(Apost); cax2=caxis;%, [0 0.4])
% B
subplot(3,2,3);
imagesc(Bpre); cax3=caxis;%, [0 0.4])
subplot(3,2,4);
imagesc(Bpost); cax4=caxis;%, [0 0.4])
% X
subplot(3,2,5);
imagesc(Xpre); cax5=caxis;%, [0 0.4])
subplot(3,2,6);
imagesc(Xpost); cax6=caxis;%, [0 0.4])


subplot(3,2,1);caxis([min(cax1(1),cax2(1)),max(cax1(2),cax2(2))]); colorbar; title('Apre')
subplot(3,2,2);caxis([min(cax1(1),cax2(1)),max(cax1(2),cax2(2))]); colorbar; title('Apost')

subplot(3,2,3);caxis([min(cax3(1),cax4(1)),max(cax3(2),cax4(2))]); colorbar; title('Bpre')
subplot(3,2,4);caxis([min(cax3(1),cax4(1)),max(cax3(2),cax4(2))]); colorbar; title('Bpost')

subplot(3,2,5);caxis([min(cax5(1),cax6(1)),max(cax5(2),cax5(2))]); colorbar; title('Xpre')
subplot(3,2,6);caxis([min(cax5(1),cax6(1)),max(cax6(2),cax6(2))]); colorbar; title('Xpost')


%% Statistics

% load data A per subject
condname    = 'X';
fname       = 'cohspctrm';
condidx     = [1 4; 2 5; 3 6];
if condname == 'A'
    row = 1;
elseif condname == 'B'
    row = 2;
elseif condname == 'X'
    row = 3;
end
for i = 1:length(subjects)
    
    disp(strcat('Subject triangulation:', int2str(i)))
    load(fullfile('/project/3012026.13/jansch/',strcat(subjects{i}, '_coh.mat')))
    Pre            = coh{condidx(row,1)}.(fname);
    Post           = coh{condidx(row,2)}.(fname);
    triang         = tril(Pre,-1);
    dataPre(:,i)   = Pre(triang>0);
    dataPost(:,i)  = Post(triang>0);

end

% en ook zo met B X pre en post

nsubj       = length(subjects);
cfg         = [];
design      = [ones(1,nsubj) ones(1,nsubj)*2; 1:nsubj 1:nsubj];
cfg.ivar    = 1;cfg.uvar=2;

stat        = ft_statfun_depsamplesT(cfg, [dataPre dataPost], design);

Tstat       = zeros(size(Pre));
Tstat(triang>0) = stat.stat;

statname = fullfile('/project/3012026.13/jansch/', sprintf('%s_%s_%s', 'Tstats', condname, fname));
save(statname, 'Tstat', 'fname');


%% Load and plot statistics -  Take sorted_idx from above

load(fullfile('/project/3012026.13/jansch/', 'Tstats_A_cohspctrm.mat'))
Tstat       = Tstat+Tstat';
Apre_post   = Tstat(sorted_idx,sorted_idx);
load(fullfile('/project/3012026.13/jansch/', 'Tstats_B_cohspctrm.mat'))
Tstat       = Tstat+Tstat';
Bpre_post   = Tstat(sorted_idx,sorted_idx);
load(fullfile('/project/3012026.13/jansch/', 'Tstats_X_cohspctrm.mat'))
Tstat       = Tstat+Tstat';
Xpre_post   = Tstat(sorted_idx,sorted_idx);

figure;
% A pre post
subplot(3,1,1);
imagesc(Apre_post); cax1=caxis;
% B pre post
subplot(3,1,2);
imagesc(Bpre_post); cax2=caxis;
% X pre post
subplot(3,1,3);
imagesc(Bpre_post); cax3=caxis;

subplot(3,1,1);caxis([min([cax1(1),cax2(1),cax3(1)]),max([cax1(2),cax2(2),cax3(1)])]); colorbar; title('A (pre vs post)')
subplot(3,1,2);caxis([min([cax1(1),cax2(1),cax3(1)]),max([cax1(2),cax2(2),cax3(1)])]); colorbar; title('B (pre vs post)')
subplot(3,1,3);caxis([min([cax1(1),cax2(1),cax3(1)]),max([cax1(2),cax2(2),cax3(1)])]); colorbar; title('X (pre vs post)')



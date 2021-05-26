%% Pipeline to initiate RSA on source reconstructed data
%%
addpath /home/common/matlab/fieldtrip/
addpath /home/common/matlab/fieldtrip/qsub/
addpath /home/common/matlab/fieldtrip/private/
addpath /home/common/matlab/fieldtrip/scripts_source_connectivity/
addpath /project/3012026.13/scripts_RT/
ft_defaults

root_dir     = '/project/3012026.13/processed/';
save_dir     = '/project/3012026.13/processed_RT/source reconstruction/';

subjects = {'sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-022','sub-023',...
    'sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-031'...
    'sub-032','sub-033','sub-034','sub-035','sub-036','sub-037','sub-038'};

%% Call the source RSA function

% Needs filterDataStory & my_data from create_virtual_dipoles.m

load(fullfile(save_dir, 'atlas_grid.mat'), 'atlas_grid')

% Apply searchlight function
cfg                 = [];
cfg.timewidth       = 0.250; % Width of shifting timewindow
cfg.timesteps       = 0.125; % Steps of shifting timewindow
cfg.time            = my_data{1}.time;
cfg.atlas           = atlas_grid;
cfg.parcelrm        = []; % These parcels for L and R hemisphere & medial wall are removed
cfg.searchspace     = 'yes'; % Not included in rt_searchlight yet
cfg.searchtime      = 'yes'; % Not included in rt_searchlight yet
cfg.avgovertrials   = 'yes'; % Not included in rt_searchlight yet
cfg.avgovertime     = 'no'; % Not included in rt_searchlight yet
cfg.avgoverspace    = 'no'; % Not included in rt_searchlight yet
dataRSA             = rt_sourcersa(cfg, filterDataStory);
        
save(fullfile(save_dir, subjects{3}, 'dataRSA.mat'), 'dataRSA')


%% Create source structure based on RSA data and average over subjects

ftpath   = '/home/common/matlab/fieldtrip'; % this is the path to FieldTrip at Donders
load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d7point5mm'));
template_grid = sourcemodel;
clear sourcemodel;
load(fullfile(save_dir, 'atlas_grid.mat'), 'atlas_grid')

ignoreSubs = [4 10 12 20];
normidx = 1;
for normsub = 1:10
    
    disp(strcat('Calculating general filter: ', int2str(normsub)))
    disp('Loading...')
    
    if ~ismember(normsub, ignoreSubs)

        load(fullfile(save_dir, subjects{normsub}, 'dataRSA.mat'), 'dataRSA')
        load(fullfile(save_dir, subjects{normsub}, 'source_filter.mat'), 'source')
        rsaparcellation = nan(length(atlas_grid.parcellation1D),size(dataRSA.RSM, 2));
        
        disp(strcat('Calculating parcel:...'))
        for parcel = 1:length(dataRSA.parcels)

            idxparcel                    = atlas_grid.parcellation1D == dataRSA.parcels(parcel);
            parcel_matrix                = ones(size(rsaparcellation(idxparcel,:))).*dataRSA.RSM(parcel,:);
            rsaparcellation(idxparcel,:) = parcel_matrix;

        end

        source_rsa                  = source;
        source_rsa.rsaparcellation  = rsaparcellation(:,1);
        source_rsa.parcels          = atlas_grid.parcellation1D;
        source_rsa.pos              = template_grid.pos;
        source_rsa.inside           = template_grid.inside;

        source_rsaAll{normidx,1}    = source_rsa;
        
        normidx = normidx + 1;

    end
end

cfg           = [];
cfg.parameter = 'rsaparcellation';
source_rsaAVG = ft_sourcegrandaverage(cfg, source_rsaAll{:,1});

mri = ft_read_mri('/home/common/matlab/fieldtrip/template/anatomy/single_subj_T1.nii');
        
cfg                         = [];
cfg.interpmethod            = 'nearest';
cfg.parameter               = 'rsaparcellation';
source_rsaInt               = ft_sourceinterpolate(cfg, source_rsaAVG, mri);
source_rsaInt.rsaparcellation(source_rsaInt.rsaparcellation == 0) = NaN;

% plot multiple 2D axial slices
cfg = [];
cfg.method        = 'surface';
cfg.surfinflated  = 'surface_inflated_both.mat';
cfg.funparameter  = 'tissue';
%cfg.maskparameter = cfg.funparameter;
%cfg.funcolorlim   = [0 0.08];
%cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';
cfg.camlight       = 'no';
figure;ft_sourceplot(cfg, atlas);
set(gcf,'color','w')
%ft_hastoolbox('brewermap', 1);         
%colormap(brewermap(64,'OrRd'))
view([-20 30])


%% Create channel X timewindow plot

[nr,nc] = size(dataRSA.RSM);
pcolor([dataRSA.RSM nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))


%% Create line plot with channel x timewindow trial data

meanLine = mean(dataRSA.RSM,1);
timeRSM = [1:1:15];

figure('Renderer', 'painters', 'Position', [10 10 500 200])
hold on
plot(timeRSM, dataRSA.RSM, 'color', [0.5 0.5 0.5]);
plot(timeRSM, meanLine, 'color', [0.6350 0.0780 0.1840]);
%ylim([-0.04 0.04])
xlim([1 15])
xlabel('Time-window')
ylabel('Similarity (r)')
title('Similarity for each channel over time')
set(gcf,'color','w')


%% Test atlas

% This was used to see where a selection of parcels will show up in source
% space. source_rsaInt can then be used for plotting, including the
% respective mask

atlas_roi = {
'R_18_B05_01'     
'R_18_B05_02'     
'R_18_B05_03'     
'R_18_B05_04'     
'R_18_B05_05'     
'R_18_B05_06'     
'R_18_B05_07'     
'L_18_B05_01'     
'L_18_B05_02'     
'L_18_B05_03'     
'L_18_B05_04'     
'L_18_B05_05'     
'L_18_B05_06'     
'L_18_B05_07'}

% Check which brodmann areas carries which label
idxRois = [];
for entries = 1:length(source_rsaInt.parcellationlabel)
    for rois = 1:length(atlas_roi)
        if strcmp(source_rsaInt.parcellationlabel{entries}, atlas_roi{rois})
            idxRois = [idxRois; entries];
        end
    end 
end

allidxs = [];
for myrois = 1:length(idxRois)
    parcidx = find(source_rsaInt.parcellation == idxRois(myrois));
    allidxs = [allidxs; parcidx];
end

parcellation2 = source_rsaInt.parcellation;
parcellation2(allidxs) = 500;
parcellation2(parcellation2 < 500) = 0;
mymask = parcellation2 == 500;
source_rsaInt.mymask = mymask;



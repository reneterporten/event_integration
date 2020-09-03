%% Create virtual dipoles for each parcel

%% Default path

addpath /home/common/matlab/fieldtrip/
addpath /home/common/matlab/fieldtrip/qsub/
addpath /project/3012026.13/scripts_RT/
ft_defaults

root_dir     = '/project/3012026.13/processed/';
save_dir     = '/project/3012026.13/processed_RT/source reconstruction/';

subjects = {'sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-022','sub-023',...
    'sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-031'...
    'sub-032','sub-033','sub-034','sub-035','sub-036','sub-037','sub-038'};


%% Source analysis - Create spatial filters based on the whole data

load(fullfile(save_dir, subjects{2}, 'headmodel.mat'), 'headmodel')
load(fullfile(save_dir, subjects{2}, 'grid.mat'), 'grid')

% Load Data
tlck = rt_mytimelock_cov(root_dir, subjects{2});

[u,s,v] = svd(tlck.cov);
d       = -diff(log10(diag(s)));
d       = d./std(d);
kappa   = find(d>5,1,'first');

figure;
semilogy(diag(s),'o-');

cfg             = [];
cfg.grid        = grid;
cfg.headmodel   = headmodel;
cfg.channel     = {'MEG'};
cfg.grad        = tlck.grad;
sourcemodel_lf  = ft_prepare_leadfield(cfg, tlck);

 % Engange the source analysis for specified ROIs only
cfg                     = [];
cfg.method              = 'lcmv';
cfg.grad                = tlck.grad;
cfg.lcmv.kappa          = kappa;
cfg.lcmv.keepfilter     = 'yes';
cfg.lcmv.fixedori       = 'yes';
cfg.lcmv.weightnorm     = 'unitnoisegain';
cfg.lcmv.lambda         = '5%';
cfg.lcmv.kurtosis       = 'yes';
cfg.headmodel           = headmodel;
cfg.sourcemodel         = sourcemodel_lf;
source                  = ft_sourceanalysis(cfg, tlck);


%% Plotting the kurtosis (the spikeness of the data)

load(fullfile(save_dir, subjects{2}, 'alignedmri.mat'), 'mri')

cfg                 = [];
cfg.parameter       = 'kurtosis';
source_interp       = ft_sourceinterpolate(cfg, testSource, mri);

cfg                 = [];
cfg.funparameter    = 'kurtosis';
cfg.method          = 'slice'; % orthogonal slices with crosshairs at peak (default anyway if not specified)
figure;ft_sourceplot(cfg, source_interp);

clear source_interp


%% Apply general filter to sensor data

% Load Data
subj                = subjects{2};
% Apply timelock analysis
my_data             = rt_mytimelock_cov2(root_dir, subj);

filterDataStory = cell(length(my_data)*length(my_data{1}.trialinfo),1);
running_idx     = 1;
for stor = 1:length(my_data)
    
    disp(strcat('Story:', int2str(stor)))
    
    my_data_sel         = my_data{stor};

    % First create a data structure that can store all filter data per
    % location x channel coefficients. 
    filterData  = zeros(size(source.pos, 1), numel(my_data_sel.label));

    % Subsequently take the filter data and store them in the pre-allocated
    % variable.
    filterData(source.inside,:) = cat(1,source.avg.filter{:});
    
    for trials = 1:length(my_data_sel.trialinfo)
        % In order to obtain per location the dipole level fourier data, I multiply
        % the sensor timelocked data with the spatial filter obtained from the source
        % analysis.
        filterDataTrial                 = filterData*squeeze(my_data_sel.trial(trials,:,:));
        filterDataStory{running_idx}    = filterDataTrial;
        running_idx                     = running_idx + 1;
    end
    
end


%% Prepare atlas alignment

load(fullfile(save_dir, subjects{1}, 'grid.mat'), 'grid')

atlas = ft_read_atlas('/project/3012026.13/scripts_RT/atlas_subparc374avg_volumetric.mat');
atlas = ft_convert_units(atlas, 'mm');

% The parcellation of the volumetric atlas gets interpolated onto the
% subject specific sourcemodel (grid)
cfg                         = [];
cfg.interpmethod            = 'nearest';
cfg.parameter               = 'parcellation';
atlas_grid                  = ft_sourceinterpolate(cfg, atlas, grid);
atlas_grid.parcellation1D   = atlas_grid.parcellation(:);

save(fullfile(save_dir, subjects{1}, 'atlas_grid.mat'), 'atlas_grid')


%% Calculate number of grid points for each parcel

% Depending on the template_grid that was used for the subject specific
% sourcemodel, the amount of gird points might change. This has an
% influence on how many grid points are associated with individual parcels
% of the atlas. This routine plots the distribution of grid points for each
% parcel

uniqueParc = unique(atlas_grid.parcellation1D);
uniqueParc = uniqueParc(2:end);

allNumPar = zeros(length(uniqueParc),1);
for par = 1:length(uniqueParc)
    
    currentPar = uniqueParc(par);
    numPar = sum(atlas_grid.parcellation1D == currentPar);
    allNumPar(par) = numPar;
    
end

plot(uniqueParc, allNumPar, '.')

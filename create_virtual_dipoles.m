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


%% Estimate covariance

% Estimate covariance per subject
for covsub = 1:10
    qsubfeval('rt_mytimelock_cov', root_dir, subjects{covsub}, 'memreq', 12*1024^3, 'timreq', 24*60, 'batchid', subjects{covsub});
end


%% Source analysis - Create spatial filters based on the whole data

for sourcesub = 1:10
    
    disp(strcat('Calculating general filter: ', int2str(sourcesub)))
    disp('Loading...')
  
    load(fullfile(save_dir, subjects{sourcesub}, 'headmodel.mat'), 'headmodel')
    load(fullfile(save_dir, subjects{sourcesub}, 'grid.mat'), 'grid')
    load(fullfile(save_dir, subjects{sourcesub}, 'tlck_filter.mat'), 'tlck')

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
    cfg.normalize   = 'yes';
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
    
    save(fullfile(save_dir, subjects{sourcesub}, 'source_filter.mat'), 'source')
    clear source sourcemodel_lf tlck headmodel grid

end


%% Plotting the kurtosis (the spikeness of the data)

load(fullfile(save_dir, subjects{1}, 'alignedmri.mat'), 'mri')

source.pos = template_grid.pos;
source.inside = template_grid.inside;

cfg                 = [];
cfg.parameter       = 'pow';
source_interp       = ft_sourceinterpolate(cfg, source, mri);

cfg                 = [];
cfg.funparameter    = 'pow';
cfg.method          = 'slice'; % orthogonal slices with crosshairs at peak (default anyway if not specified)
figure;ft_sourceplot(cfg, source_interp);

clear source_interp


%% Calculate timelocked data for all stories

for tlcksub = 1:10
    qsubfeval('rt_mytimelock_cov2', root_dir, subjects{tlcksub}, 'memreq', 40*1024^3, 'timreq', 40*60, 'batchid', subjects{tlcksub});
end


%% Apply general filter to sensor data

ignoreSubs = [4 10 12 20];
for virtsub = 4:10
    
    if ~ismember(virtsub, ignoreSubs)
        
        disp(strcat('Calculating general filter: ', int2str(virtsub)))
        disp('Loading...')

        load(fullfile(save_dir, subjects{virtsub}, 'tlck_structure.mat'), 'my_data')
        load(fullfile(save_dir, subjects{virtsub}, 'source_filter.mat'), 'source')

        filterDataStory = cell(length(my_data)*length(my_data{1}.trialinfo),1);
        running_idx     = 1;
        for stor = 1:length(my_data)

            disp(strcat('Cond.:', int2str(stor)))

            my_data_sel         = my_data{stor};

            % First create a data structure that can store all filter data per
            % location x channel coefficients. 
            filterData  = zeros(size(source.pos, 1), numel(my_data_sel.label));

            % Subsequently take the filter data and store them in the pre-allocated
            % variable.
            filterData(source.inside,:) = cat(1,source.avg.filter{:});

            for trials = 1:size(my_data_sel.trialinfo,1)
                % In order to obtain per location the dipole level fourier data, I multiply
                % the sensor timelocked data with the spatial filter obtained from the source
                % analysis.
                filterDataTrial                 = filterData*squeeze(my_data_sel.trial(trials,:,:));
                filterDataStory{running_idx}    = filterDataTrial;
                running_idx                     = running_idx + 1;
            end

        end

        % Call the source RSA function
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

        save(fullfile('/project/3012026.13/processed_RT/source reconstruction/', subjects{virtsub}, 'dataRSA.mat'), 'dataRSA')

        clear filterDataStory my_data_sel my_data filterData filterDataTrial atlas_grid cfg dataRSA
    
    end

end


%% Prepare atlas alignment

ftpath   = '/home/common/matlab/fieldtrip'; % this is the path to FieldTrip at Donders
load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d7point5mm'));
%load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d4mm'));
template_grid = sourcemodel;
clear sourcemodel;

template_grid = ft_convert_units(template_grid, 'mm');

atlas = ft_read_atlas('/home/common/matlab/fieldtrip/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');
atlas = ft_convert_coordsys(atlas, 'ctf');

mri = ft_read_mri('/home/common/matlab/fieldtrip/template/anatomy/single_subj_T1.nii');
        
cfg                         = [];
cfg.interpmethod            = 'nearest';
cfg.parameter               = 'tissue';
source_rsaInt               = ft_sourceinterpolate(cfg, atlas, mri);

% The parcellation of the volumetric atlas gets interpolated onto the
% subject specific sourcemodel (grid)
cfg                         = [];
cfg.interpmethod            = 'nearest';
cfg.parameter               = 'tissue';
atlas_grid                  = ft_sourceinterpolate(cfg, atlas, template_grid);
atlas_grid.parcellation1D   = atlas_grid.tissue(:);

save(fullfile(save_dir, 'atlas_grid.mat'), 'atlas_grid')


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

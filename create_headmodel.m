%% Create headmodel for source RSA

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

%% MRI data location

for subseg = 1:10 %length(subjects)
    
    dicoms = fullfile('/project/3012026.13/anatomy/', subjects{subseg}, '/preproc/', strcat(subjects{subseg}, '_mri.mgz'));

    mri = ft_read_mri(dicoms);

    cfg         = [];
    cfg.method  = 'flip';
    mri         = ft_volumereslice(cfg, mri);

    %Define coordinates
    cfg             = [];
    cfg.method      = 'interactive';
    cfg.coordsys    = 'spm';
    cfg.parameter   = 'anatomy';
    [mri]           = ft_volumerealign(cfg, mri);

    cfg             = [];
    cfg.resolution  = 1;
    cfg.dim         = [256 256 256];
    mri             = ft_volumereslice(cfg, mri);

    cfg             = [];
    cfg.method      = 'interactive';
    cfg.coordsys    = 'ctf';
    mri             = ft_volumerealign(cfg, mri);

    cfg             = [];
    seg             = ft_volumesegment(cfg, mri);

    % make a figure of the mri and segmented volumes to visually inspect the
    % segmentation
    %segmentedmri           = seg;
    %segmentedmri.transform = mri.transform;
    %segmentedmri.anatomy   = mri.anatomy;
    %cfg                    = [];
    %cfg.funparameter       = 'gray';
    %ft_sourceplot(cfg, segmentedmri);

    mkdir(fullfile(save_dir, subjects{subseg}))
    save(fullfile(save_dir, subjects{subseg}, 'segmentedmri.mat'), 'seg')
    save(fullfile(save_dir, subjects{subseg}, 'alignedmri.mat'), 'mri')

end


%% Create headmodel

for headsub = 1:10
    
    load(fullfile(save_dir, subjects{headsub}, 'segmentedmri.mat'), 'seg')
    
    % compute the subject's headmodel/volume conductor model
    cfg                = [];
    cfg.method         = 'singleshell';
    headmodel          = ft_prepare_headmodel(cfg, seg);

    save(fullfile(save_dir, subjects{headsub}, 'headmodel.mat'), 'headmodel')

end


%% Create leadfield based on template grid

ftpath   = '/home/common/matlab/fieldtrip'; % this is the path to FieldTrip at Donders
load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d7point5mm'));
template_grid = sourcemodel;
clear sourcemodel;

template_grid = ft_convert_units(template_grid, 'mm');

for gridsub = 1:10

    load(fullfile(save_dir, subjects{gridsub}, 'alignedmri.mat'), 'mri')
    
    % create the subject specific grid, using the template grid that has just been created
    cfg           = [];
    cfg.warpmni   = 'yes';
    cfg.template  = template_grid;
    cfg.nonlinear = 'yes';
    cfg.mri       = mri;
    cfg.unit      ='mm';
    grid          = ft_prepare_sourcemodel(cfg);

    save(fullfile(save_dir, subjects{gridsub}, 'grid.mat'), 'grid')

end


function [headmodel, sourcemodel] = rt_headsource(subj, varargin)

% Function to define the headmodel

if nargin<1 || isempty(subj)
  subj = 'sub-004';
end

mridata     = ft_getopt(varargin, 'cfgpreproc', fullfile('/project/3012026.13/jansch/', strcat(subj, '_mrialigned.mat')));
template    = ft_getopt(varargin, 'template', '/home/common/matlab/fieldtrip/template/sourcemodel/standard_sourcemodel3d5mm');
saveflag    = ft_getopt(varargin, 'saveflag', true);
savepath    = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
savename    = ft_getopt(varargin, 'savename', 'headsource'); 


%% Load and segment MRI

load(mridata)

cfg             = [];
seg             = ft_volumesegment(cfg, mri);


%% Create sourcemodel

load(template)
template_grid = sourcemodel;
template_grid = ft_convert_units(template_grid, 'mm');
clear sourcemodel

% ft_prepare_sourcemodel
cfg           = [];
cfg.warpmni   = 'yes';
cfg.template  = template_grid;
cfg.nonlinear = 'yes';
cfg.mri       = mri;
cfg.unit      ='mm';
sourcemodel   = ft_prepare_sourcemodel(cfg);


%% compute the subject's headmodel/volume conductor model

cfg             = [];
cfg.method      = 'singleshell';
headmodel       = ft_prepare_headmodel(cfg, seg);


%% Save data

if saveflag
    fname = fullfile(savepath, sprintf('%s_%s', subj, savename));
    save(fname, 'sourcemodel', 'headmodel');
end



function [headmodel, sourcemodel] = rt_headsource(cfg, varargin)

% Function to define the headmodel

if nargin<1 || isempty(subj)
  subj = 'sub-004';
end

mridata     = ft_getopt(varargin, 'cfgpreproc', fullfile('/project/3012026.13/jansch/', strcat(subj, '_mri.mat')));
saveflag    = ft_getopt(varargin, 'saveflag', true);
savepath    = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
savename    = ft_getopt(varargin, 'savename', 'headsource'); 


%% Load and segment MRI

load(mridata)

cfg             = [];
seg             = ft_volumesegment(cfg, mri);


%% Create sourcemodel

% ft_prepare_sourcemodel


%% compute the subject's headmodel/volume conductor model

cfg             = [];
cfg.method      = 'singleshell';
headmodel       = ft_prepare_headmodel(cfg, seg);


%% Save data

if saveflag
    fname = fullfile(savepath, sprintf('%s_%s', subj, savename));
    save(fname, 'sourcemodel', 'headmodel');
end



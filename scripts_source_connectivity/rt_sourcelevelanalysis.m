function rt_sourcelevelanalysis(cfg, varargin)

% Function to define leadfield and estimate activity at source level
% Returns connectivity measures between the hippocampus and the mPFC

if nargin<1 || isempty(subj)
  subj = 'sub-004';
end

saveflag    = ft_getopt(varargin, 'saveflag', false);
sensdata    = ft_getopt(varargin, 'sensdata', fullfile('/project/3012026.13/jansch/', strcat(subj, '_nocomb.mat')));
headdata    = ft_getopt(varargin, 'sensdata', fullfile('/project/3012026.13/jansch/', strcat(subj, '_headsource.mat')));


%% Load time-frequency data (output from rt_sensorlevelanalysis) & headmodel

load(sensdata) % 1x6 freq structure
load(headdata) % Loads headmodel and segmented mri


%% Prepare leadfield

cfg                     = [];
cfg.grad                = freq(1).grad;
cfg.headmodel           = headmodel;
cfg.channel             = {'MEG'};
sourcemodel             = ft_prepare_leadfield(cfg);


%% Source analysis

cfg                     = [];
cfg.method              = 'dics';
cfg.frequency           = 18;
cfg.sourcemodel         = sourcemodel;
cfg.headmodel           = headmodel;
cfg.dics.projectnoise   = 'yes';
cfg.dics.lambda         = 0;

sourcePost_nocon        = ft_sourceanalysis(cfg, freq(1));
disp('done')









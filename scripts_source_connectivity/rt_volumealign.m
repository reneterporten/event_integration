function [mri] = rt_volumealign(subj, varargin)

% Function to align the structural MRI for further analyses

if nargin<1 || isempty(subj)
  subj = 'sub-004';
end

mripath     = ft_getopt(varargin, 'cfgpreproc', fullfile('/project/3012026.13/anatomy/', subj, '/preproc/', strcat(subj, '_mri.mgz')));
saveflag    = ft_getopt(varargin, 'saveflag', true);
savepath    = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
savename    = ft_getopt(varargin, 'savename', 'mrialigned'); 


%% Align and segment MRI

mri = ft_read_mri(mripath);

cfg             = [];
cfg.dim         = [256 256 256];
mri             = ft_volumereslice(cfg, mri);

%Define coordinates
cfg             = [];
cfg.method      = 'interactive';
cfg.coordsys    = 'ctf';
cfg.parameter   = 'anatomy';
[mri]           = ft_volumerealign(cfg, mri);


%% Save data

if saveflag
    fname = fullfile(savepath, sprintf('%s_%s', subj, savename));
    save(fname, 'mri');
end



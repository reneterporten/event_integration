function [atlas_grid] = rt_atlasalign(varargin)

% Function to align the brainnetome atlas with a tempalte grid

brainatlas  = ft_getopt(varargin, 'template', '/home/common/matlab/fieldtrip/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');
template    = ft_getopt(varargin, 'template', '/home/common/matlab/fieldtrip/template/sourcemodel/standard_sourcemodel3d5mm');
saveflag    = ft_getopt(varargin, 'saveflag', true);
savepath    = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
savename    = ft_getopt(varargin, 'savename', 'brainnetome_atlas_grid'); 


%% Load atlas and interpolate atlas onto template grid

% read the atlas
atlas               = ft_read_atlas(brainatlas);

% load the template sourcemodel
load(template)
template_grid       = sourcemodel;
template_grid       = ft_convert_units(template_grid, 'mm');
clear sourcemodel

% Source-interpolate atlas onto tempalte_grid
cfg                 = [];
cfg.interpmethod    = 'nearest';
cfg.parameter       = 'tissue';
cfg.unit            = 'mm';
atlas_grid          = ft_sourceinterpolate(cfg, atlas, template_grid);


%% Save aligned atlas

if saveflag
    fname = fullfile(savepath, sprintf('%s', savename));
    save(fname, 'atlas_grid');
end



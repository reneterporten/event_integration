function rt_visualizeroi(varargin)

% Function for plotting atlas specific parcels onto template mri

brainatlas      = ft_getopt(varargin, 'brainatlas', '/home/common/matlab/fieldtrip/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');
template_mri    = ft_getopt(varargin, 'brainatlas', '/project/3012026.13/scripts_RT/git_event_integration/fieldtrip/template/anatomy/single_subj_T1_1mm.nii');
atlasrois       = ft_getopt(varargin, 'atlasrois', 'all'); % Either all or cell array with ROIs
plotmethod      = ft_getopt(varargin, 'plotmethod', 'ortho');
plotcolors      = ft_getopt(varargin, 'plotcolors', 'Paired');


%% Load atlas and plot

orig_atlas  = ft_read_atlas(brainatlas);
mri         = ft_read_mri(template_mri);

cfg             = [];
cfg.coordsys    = 'mni';
cfg.dim         = mri.dim;
mri             = ft_volumereslice(cfg, mri);

cfg                 = [];
cfg.parameter       = 'tissue';
cfg.interpmethod    = 'nearest';
atlas_int           = ft_sourceinterpolate(cfg, orig_atlas , mri);
atlas_int.coordsys  = 'mni';

% Select ROIs
if strcmp(atlasrois, 'all')
    idx = 1:numel(orig_atlas.tissuelabel);
else
    idx = [];
    for lab = 1:numel(orig_atlas.tissuelabel)
        for roirun = 1:numel(atlasrois)
            if contains(orig_atlas.tissuelabel{lab}, atlasrois(roirun))
                idx = [idx; lab];
            end
        end
    end
end

% Sourceplot ROI of atlas onto tempalte MRI
ft_hastoolbox('brewermap', 1); 
cfg                 = [];
cfg.method          = plotmethod;
cfg.funparameter    = 'tissue';
cfg.funcolormap     = colormap(flipud(brewermap(numel(unique(orig_atlas.tissue)), plotcolors)));
cfg.atlas           = orig_atlas;
cfg.roi             = orig_atlas.tissuelabel(idx);
ft_sourceplot(cfg, atlas_int)
% Only orthogonal presentation allows for funcolormap, for all other cases:
if ~strcmp(plotmethod, 'ortho')
    colormap(flipud(brewermap(numel(unique(orig_atlas.tissue)), plotcolors)))
end



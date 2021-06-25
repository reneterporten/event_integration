function rt_visualizeroi(varargin)

% Function for plotting atlas specific parcels onto template brain surface

atlasgrid  = ft_getopt(varargin, 'atlasgrid', '/project/3012026.13/jansch/brainnetome_atlas_grid.mat');
brainatlas  = ft_getopt(varargin, 'brainatlas', '/home/common/matlab/fieldtrip/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');


%% Load atlas and plot

load(atlasgrid)

atlaslabel      = cell(numel(atlas_grid.tissue),1);
origatlasidx    = reshape(atlas_grid.tissue, [numel(atlas_grid.tissue),1]);
for labs = 1:numel(origatlasidx)
    if origatlasidx(labs) == 0
        atlaslabel{labs,1} = '';
    else
        atlaslabel{labs,1} = atlas_grid.tissuelabel{origatlasidx(labs,1)};
    end
end

atlaslabel = reshape(atlaslabel, [size(atlas_grid.tissue)]);


cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas_grid)






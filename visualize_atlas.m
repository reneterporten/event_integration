% Visualize structures of the atlas

addpath /home/common/matlab/fieldtrip/
addpath /home/common/matlab/fieldtrip/qsub/
addpath /project/3012026.13/scripts_RT/
ft_defaults

%%
% Broadmann areas for mpfc
% 12 Orbitofrontal area (used to be part of BA11, refers to the area between the superior frontal gyrus and the inferior rostral sulcus)
% 25 Subgenual area (part of the Ventromedial prefrontal cortex)
% 32 Dorsal anterior cingulate cortex
% 33 Part of anterior cingulate cortex
% 24 Ventral anterior cingulate cortex

brainnetome = ft_read_atlas('/home/common/matlab/fieldtrip/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');

imagesc(brainnetome.tissue(:,:,68))

idx_roi = [];
for labels = 1:length(brainnetome.tissuelabel)
    
    disp(brainnetome.tissuelabel{labels})
    
    %if contains(brainnetome.tissuelabel{labels}, 'Hipp')
    %    disp(brainnetome.tissuelabel{labels})
    %    idx_roi =[idx_roi, labels];
    %end
    
    
end

my_brainnetome = brainnetome;
my_brainnetome.tissue = my_brainnetome.tissue(:);

not_roi = find(~ismember(my_brainnetome.tissue(:), idx_roi));
is_roi = find(ismember(my_brainnetome.tissue(:), idx_roi));
my_brainnetome.tissue(not_roi) = 0;
my_brainnetome.tissue(is_roi) = 10;




mri = ft_read_mri('/home/common/matlab/fieldtrip/template/anatomy/single_subj_T1_1mm.nii');


cfg                 = [];
cfg.parameter       = 'tissue';
cfg.interpmethod    = 'nearest';
source_interp       = ft_sourceinterpolate(cfg, brainnetome, mri);
source_interp.tissue(source_interp.tissue == 0) = NaN;


cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'tissue';
cfg.maskparameter = cfg.funparameter;
%cfg.funcolorlim   = [0 0.08];
%cfg.opacitylim    = [0.0 1.2];
%cfg.opacitymap    = 'rampup';
cfg.camlight       = 'no';
figure;ft_sourceplot(cfg, source_interp);

set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(brewermap(246,'Paired'))




function rt_collectsourcedata(suff, varargin)

% Function that calculates statistics based on connectivity data

saveflag    = ft_getopt(varargin, 'saveflag', false);
redcomp     = ft_getopt(varargin, 'redcomp', true);
savepath    = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
datadir     = ft_getopt(varargin, 'datadir', '/project/3012026.13/jansch/');
savename    = ft_getopt(varargin, 'savename', 'coherence'); 
suff        = ft_getopt(varargin, 'suff', '_coh.mat');
atlasgrid   = ft_getopt(varargin, 'headdata', fullfile('/project/3012026.13/jansch/', 'brainnetome_atlas_grid.mat'));

cd(datadir);
load(atlasgrid)

cfg1 = [];
cfg1.appenddim = 'rpt';

cfg2 = [];
cfg2.avgoverrpt = 'yes';

d = dir(sprintf('sub*%s*',suff));
for k = 1:numel(d)
  
  disp(strcat('Subject:', int2str(k)))
    
  data = load(d(k).name);
  
  istimelock = ft_datatype(data.cohpre, 'timelock');
  isfreq     = ft_datatype(data.cohpre, 'freq');
  
  varname = fieldnames(data);
  for m = 1:numel(varname)
    disp(strcat('Data pre or post:', int2str(m)))
    % In case coherence of components needs to be reduced to 1
    % !!! Should include routine to check actual number of components in
    % !!! the data
    if redcomp
        % Component reduction from 5 to 1
        red_size        = size(data.(varname{m}).cohspctrm)/5;
        red_cohspctrm   = zeros(red_size);
        x_red           = 1;
        y_red           = 1;
        for y = 1:5:length(data.(varname{m}).cohspctrm)
            for x = 1:5:length(data.(varname{m}).cohspctrm)
                msnip       = data.(varname{m}).cohspctrm(x:(x+4), y:(y+4));
                el_msnip    = numel(msnip);
                msnip       = sum(sum(msnip))/el_msnip;
                red_cohspctrm(x_red,y_red) = msnip;
                x_red = x_red +1;
            end
            y_red = y_red +1;
            x_red = 1;
        end
        data.(varname{m}).cohspctrm = red_cohspctrm;
        data.(varname{m}).label     = atlas_grid.tissuelabel';
    end
    F{m,k} = data.(varname{m});
  end
  clear data
end

if istimelock
  for m = 1:size(F,1)
    Fcoh{m} = ft_selectdata(cfg2, ft_appendtimelock(cfg1, F{:,m}));
  end
elseif isfreq
   for m = 1:size(F,1) 
     Fcoh{m} = ft_selectdata(cfg2, ft_appendfreq(cfg1, F{:,m}));
   end
end

clear F

subj = {d.name}';
for k = 1:numel(subj)
  subj{k} = subj{k}(1:7);
end


%% Save variables

if saveflag
    fname = fullfile(savepath, sprintf('%s_%s', 'groupdata', savename));
    save(fname, 'Fcoh');
end


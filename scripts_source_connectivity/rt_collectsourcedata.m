function rt_collectsourcedata(varargin)

% Function that calculates statistics based on connectivity data

saveflag        = ft_getopt(varargin, 'saveflag', false);
savepath        = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
datadir         = ft_getopt(varargin, 'datadir', '/project/3012026.13/jansch/');
savename        = ft_getopt(varargin, 'savename', 'coherence'); 
suff            = ft_getopt(varargin, 'suff', '_coh.mat');
connectivity    = ft_getopt(varargin, 'connectivity', 'coh');
atlasgrid       = ft_getopt(varargin, 'headdata', fullfile('/project/3012026.13/jansch/', 'brainnetome_atlas_grid.mat'));
method          = ft_getopt(varargin, 'method', 'avg'); % can also be 'stat'

cd(datadir);
load(atlasgrid)

cfg1 = [];
cfg1.appenddim = 'rpt';

cfg2 = [];
cfg2.avgoverrpt = 'yes';

d = dir(sprintf('sub*%s*',suff));
for k = 1:numel(d)
  
  disp(strcat('Subject aggregation:', int2str(k)))
    
  data_all = load(d(k).name);
  switch connectivity
      case 'coh'
          data = data_all.coh;
      case 'imcoh'
          data = data_all.imcoh;
      case 'mim'
          data = data_all.mim;
  end
  conlabel = data_all.conlabel;
  clear data_all
  
  istimelock = ft_datatype(data{1}, 'timelock');
  isfreq     = ft_datatype(data{1}, 'freq');
  
  % Loop trough each condition
  for m = 1:numel(data)
    F{m,k} = data{m};
  end
  clear data
end

if istimelock
    switch method
        case 'avg'
          for m = 1:size(F,1)
            Fcon{m} = ft_selectdata(cfg2, ft_appendtimelock(cfg1, F{:,m}));
          end
        case 'stat'
            % Do statistics
    end
elseif isfreq
    switch method
        case 'avg'
           for m = 1:size(F,1) 
             Fcon{m} = ft_selectdata(cfg2, ft_appendfreq(cfg1, F{:,m}));
           end
        case 'stat'
            % Do statistics
    end
end

clear F

subj = {d.name}';
for k = 1:numel(subj)
  subj{k} = subj{k}(1:7);
end


%% Save variables

if saveflag
    fname = fullfile(savepath, sprintf('%s_%s_%s', 'groupdata', connectivity, savename));
    save(fname, 'Fcon', 'subj', 'conlabel', 'connectivity');
end


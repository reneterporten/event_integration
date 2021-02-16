function rt_collect_FREQ(suff)

datadir = '/project/3012026.13/jansch';
cd(datadir);

cfg1 = [];
cfg1.appenddim = 'rpt';

cfg2 = [];
cfg2.avgoverrpt = 'yes';

d = dir(sprintf('sub*%s*',suff));
for k = 1:numel(d)
  k
  load(d(k).name);
  for m = 1:size(freq,1)
    % these are the indivdual pairwise correlations
 
    F{m,k} = ft_selectdata(cfg2, ft_appendfreq(cfg1, freq{m,:}));
  end
end

% ensure same time axis (subsampling misalignment may occur due to
% resampling step
for k = 1:numel(F)
  F{k}.time = F{1}.time;
  F{k}.freq = F{1}.freq;
end

subj = {d.name}';
for k = 1:numel(subj)
  subj{k} = subj{k}(1:7);
end

save(sprintf('groupdata_%s',suff), 'F', 'cmb', 'subj');


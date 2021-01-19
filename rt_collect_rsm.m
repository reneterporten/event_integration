function rt_collect_rsm(suff)

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
  for m = 1:size(tlck,1)
    % these are the indivdual pairwise correlations
 
    T{m,k} = ft_selectdata(cfg2, ft_appendtimelock(cfg1, tlck{m,:}));
  end
end

% ensure same time axis (subsampling misalignment may occur due to
% resampling step
for k = 1:numel(T)
  T{k}.time = T{1}.time;
end

subj = {d.name}';
for k = 1:numel(subj)
  subj{k} = subj{k}(1:7);
end

save(sprintf('groupdata_%s',suff), 'T', 'cmb', 'subj');


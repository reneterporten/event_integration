function stat = rt_collectdata(suff)

datadir = '/project/3012026.13/jansch';
cd(datadir);

cfg1 = [];
cfg1.appenddim = 'rpt';

cfg2 = [];
cfg2.avgoverrpt = 'yes';

d = dir(sprintf('sub*%s*',suff));
for k = 1:numel(d)
  
  tmp = load(d(k).name);
  varname = fieldnames(tmp);
  data = tmp.(varname{1});
  for m = 1:numel(data)
    % these are the indivdual pairwise correlations
    F{m,k} = data(m);
  end
end

if ft_datatype(data(1), 'timelock')
  for m = 1:size(F,1)
    Fx{m} = ft_selectdata(cfg2, ft_appendtimelock(cfg1, F{:,m}));
  end
elseif ft_datatype(data(1), 'freq')
   for m = 1:size(F,1)
     Fx{m} = ft_selectdata(cfg2, ft_appendfreq(cfg1, F{:,m}));
   end
end



subj = {d.name}';
for k = 1:numel(subj)
  subj{k} = subj{k}(1:7);
end

n = size(F,2);
design = [ones(1,n) ones(1,n)*2; 1:n 1:n];

for k = 1:numel(F)
  F{k}.time = F{1}.time;
end

cfg           = [];
cfg.design    = design;
cfg.statistic = 'depsamplesT';
cfg.numrandomization = 1000;
cfg.method = 'montecarlo';
cfg.ivar = 1;
cfg.uvar = 2;

% Apost - Apre
stat{1} = ft_timelockstatistics(cfg, F{4,:}, F{1,:});
% Bpost - Bpre
stat{2} = ft_timelockstatistics(cfg, F{5,:}, F{2,:});
% Xpost - Xpre
stat{3} = ft_timelockstatistics(cfg, F{6,:}, F{3,:});

% interaction (Bpost-Apost) - (Bpre-Apre)
for k = 1:size(F,2)
  Fdiff{1,k} = F{1,k};
  Fdiff{2,k} = F{2,k};
  Fdiff{1,k}.avg = F{5,k}.avg - F{4,k}.avg;
  Fdiff{2,k}.avg = F{2,k}.avg - F{1,k}.avg;
end
stat{4} = ft_timelockstatistics(cfg, Fdiff{1,:}, Fdiff{2,:});

save(sprintf('groupdata_%s',suff), 'stat', 'subj');


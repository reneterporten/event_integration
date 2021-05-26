function stat = rt_collectdata(suff, freqband)

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
  
  istimelock = ft_datatype(data(1), 'timelock');
  isfreq     = ft_datatype(data(1), 'freq');
  if isfreq
    for m = 1:numel(data)
      tmpcfg = [];
      tmpcfg.operation = 'log10';
      tmpcfg.parameter = 'powspctrm';
      data_(m) = ft_math(tmpcfg, data(m));
    end
    data = data_;
    clear data_;
  end
  
  for m = 1:numel(data)
    
    % these are the indivdual pairwise correlations
    F{m,k} = data(m);
  end
end

if istimelock
  for m = 1:size(F,1)
    Fx{m} = ft_selectdata(cfg2, ft_appendtimelock(cfg1, F{:,m}));
  end
elseif isfreq
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

cfg                    	= [];
cfg.method          	= 'template'; 
cfg.template         	= 'CTF275_neighb.mat';               
cfg.layout           	= 'CTF275_helmet.mat';                     
cfg.feedback         	= 'no';                            
neighbours          	= ft_prepare_neighbours(cfg, F{1}); 

cfg                   	= [];
cfg.channel           	= {'all'};
cfg.neighbours        	= neighbours;
cfg.correctm         	= 'cluster';
cfg.clusteralpha     	= 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 2;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.025;
cfg.correcttail         = 'alpha';
cfg.design              = design;
cfg.statistic           = 'depsamplesT';
cfg.numrandomization    = 1000;
cfg.method              = 'montecarlo';
cfg.ivar                = 1;
cfg.uvar                = 2;

% 1 Apre
% 2 Bpre
% 3 Xpre
% 4 Apost
% 5 Bpost
% 6 Xpost
cmb = {'Apre', 'Bpre', 'Xpre', 'Apost', 'Bpost', 'Xpost'};
cmbstat = {'Apost - Apre', 'Bpost - Bpre', 'Xpost - Xpre', 'Xpre - Apre', 'Xpre - Bpre', ...
    'Xpost - Apost','Xpost - Bpost', '(Bpost-Apost) - (Bpre-Apre)', ...
    '(Xpost-Bpost) - (Xpre-Bpre)', ...
    '(Xpost-Apost) - (Xpre-Apre)', ...
    'Xpre - ABpre', 'Xpost - ABpost', '(Xpost-ABpost) - (Xpre-ABpre)'};
if istimelock
  savename = '';
  % Apost - Apre
  stat{1} = ft_timelockstatistics(cfg, F{4,:}, F{1,:});
  % Bpost - Bpre
  stat{2} = ft_timelockstatistics(cfg, F{5,:}, F{2,:});
  % Xpost - Xpre
  stat{3} = ft_timelockstatistics(cfg, F{6,:}, F{3,:});
  % Xpre - Apre
  stat{4} = ft_timelockstatistics(cfg, F{3,:}, F{1,:});
  % Xpre - Bpre
  stat{5} = ft_timelockstatistics(cfg, F{3,:}, F{2,:});
  % Xpost - Apost
  stat{6} = ft_timelockstatistics(cfg, F{6,:}, F{4,:});
  % Xpost - Bpost
  stat{7} = ft_timelockstatistics(cfg, F{6,:}, F{5,:});
  
  % interaction (Bpost-Apost) - (Bpre-Apre)
  % interaction (Xpost-Bpost) - (Xpre-Bpre)
  % interaction (Xpost-Apost) - (Xpre-Apre)
  % interaction (Xpost-ABpost) - (Xpre-ABpre)
  for k = 1:size(F,2)
    Fdiff{1,k} = F{1,k};
    Fdiff{2,k} = F{2,k};
    Fdiff{3,k} = F{3,k};
    Fdiff{4,k} = F{4,k};
    Fdiff{5,k} = F{5,k};
    Fdiff{6,k} = F{6,k};
    Fdiff{7,k} = F{6,k}; % Xpre - ABpre
    Fdiff{8,k} = F{6,k}; % Xpost - ABpost
    F_AB{1,k} = F{1,k}; % pre
    F_AB{2,k} = F{4,k}; % post
    F_AB{1,k}.avg = (F{1,k}.avg + F{2,k}.avg)./2;
    F_AB{2,k}.avg = (F{4,k}.avg + F{5,k}.avg)./2;
    % (Bpost-Apost) - (Bpre-Apre)
    Fdiff{1,k}.avg = F{5,k}.avg - F{4,k}.avg;
    Fdiff{2,k}.avg = F{2,k}.avg - F{1,k}.avg;
    % (Xpost-Bpost) - (Xpre-Bpre)
    Fdiff{3,k}.avg = F{6,k}.avg - F{5,k}.avg;
    Fdiff{4,k}.avg = F{3,k}.avg - F{2,k}.avg;
    % (Xpost-Apost) - (Xpre-Apre)
    Fdiff{5,k}.avg = F{6,k}.avg - F{4,k}.avg;
    Fdiff{6,k}.avg = F{3,k}.avg - F{1,k}.avg;
    % (Xpost-ABpost) - (Xpre-ABpre)
    Fdiff{7,k}.avg = F{6,k}.avg - F_AB{2,k}.avg;
    Fdiff{8,k}.avg = F{3,k}.avg - F_AB{1,k}.avg;
  end
  % (Bpost-Apost) - (Bpre-Apre)
  stat{8} = ft_timelockstatistics(cfg, Fdiff{1,:}, Fdiff{2,:});
  % (Xpost-Bpost) - (Xpre-Bpre)
  stat{9} = ft_timelockstatistics(cfg, Fdiff{3,:}, Fdiff{4,:});
  % (Xpost-Apost) - (Xpre-Apre)
  stat{10} = ft_timelockstatistics(cfg, Fdiff{5,:}, Fdiff{6,:});
  
  % Xpre - ABpre
  stat{11} = ft_timelockstatistics(cfg, F{3,:}, F_AB{1,:});
  % Xpost - ABpost
  stat{12} = ft_timelockstatistics(cfg, F{6,:}, F_AB{2,:});
  % (Xpost-ABpost) - (Xpre-ABpre)
  stat{13} = ft_timelockstatistics(cfg, Fdiff{7,:}, Fdiff{8,:});
  
elseif isfreq
  cfg.parameter = 'powspctrm';
  if strcmp('theta', freqband)
      cfg.frequency = [3 7];
      cfg.avgoverfreq = 'yes';
      savename = '_theta';
  elseif strcmp('alpha', freqband)
      cfg.frequency = [8 12];
      cfg.avgoverfreq = 'yes';
      savename = '_alpha';
  elseif strcmp('beta', freqband)
      cfg.frequency = [16 20];
      cfg.avgoverfreq = 'yes';
      savename = '_beta';
  else
      cfg.frequency = 'all';
      cfg.avgoverfreq = 'no';
      savename = '';
  end
  
  % Apost - Apre
  stat{1} = ft_freqstatistics(cfg, F{4,:}, F{1,:});
  % Bpost - Bpre
  stat{2} = ft_freqstatistics(cfg, F{5,:}, F{2,:});
  % Xpost - Xpre
  stat{3} = ft_freqstatistics(cfg, F{6,:}, F{3,:});
  % Xpre - Apre
  stat{4} = ft_freqstatistics(cfg, F{3,:}, F{1,:});
  % Xpre - Bpre
  stat{5} = ft_freqstatistics(cfg, F{3,:}, F{2,:});
  % Xpost - Apost
  stat{6} = ft_freqstatistics(cfg, F{6,:}, F{4,:});
  % Xpost - Bpost
  stat{7} = ft_freqstatistics(cfg, F{6,:}, F{5,:});
  
  % interaction (Bpost-Apost) - (Bpre-Apre)
  % interaction (Xpost-Bpost) - (Xpre-Bpre)
  % interaction (Xpost-Apost) - (Xpre-Apre)
  % interaction (Xpost-ABpost) - (Xpre-ABpre)
  for k = 1:size(F,2)
    Fdiff{1,k} = F{1,k};
    Fdiff{2,k} = F{2,k};
    Fdiff{3,k} = F{3,k};
    Fdiff{4,k} = F{4,k};
    Fdiff{5,k} = F{5,k};
    Fdiff{6,k} = F{6,k};
    Fdiff{7,k} = F{6,k}; % Xpre - ABpre
    Fdiff{8,k} = F{6,k}; % Xpost - ABpost
    F_AB{1,k} = F{1,k}; % pre
    F_AB{2,k} = F{4,k}; % post
    F_AB{1,k}.powspctrm = (F{1,k}.powspctrm + F{2,k}.powspctrm)./2;
    F_AB{2,k}.powspctrm = (F{4,k}.powspctrm + F{5,k}.powspctrm)./2;
    % (Bpost-Apost) - (Bpre-Apre)
    Fdiff{1,k}.powspctrm = F{5,k}.powspctrm - F{4,k}.powspctrm;
    Fdiff{2,k}.powspctrm = F{2,k}.powspctrm - F{1,k}.powspctrm;
    % (Xpost-Bpost) - (Xpre-Bpre)
    Fdiff{3,k}.powspctrm = F{6,k}.powspctrm - F{5,k}.powspctrm;
    Fdiff{4,k}.powspctrm = F{3,k}.powspctrm - F{2,k}.powspctrm;
    % (Xpost-Apost) - (Xpre-Apre)
    Fdiff{5,k}.powspctrm = F{6,k}.powspctrm - F{4,k}.powspctrm;
    Fdiff{6,k}.apowspctrm = F{3,k}.powspctrm - F{1,k}.powspctrm;
    % (Xpost-ABpost) - (Xpre-ABpre)
    Fdiff{7,k}.powspctrm = F{6,k}.powspctrm - F_AB{2,k}.powspctrm;
    Fdiff{8,k}.powspctrm = F{3,k}.powspctrm - F_AB{1,k}.powspctrm;
  end
  % (Bpost-Apost) - (Bpre-Apre)
  stat{8} = ft_freqstatistics(cfg, Fdiff{1,:}, Fdiff{2,:});
  % (Xpost-Bpost) - (Xpre-Bpre)
  stat{9} = ft_freqstatistics(cfg, Fdiff{3,:}, Fdiff{4,:});
  % (Xpost-Apost) - (Xpre-Apre)
  stat{10} = ft_freqstatistics(cfg, Fdiff{5,:}, Fdiff{6,:});
  
  % Xpre - ABpre
  stat{11} = ft_freqstatistics(cfg, F{3,:}, F_AB{1,:});
  % Xpost - ABpost
  stat{12} = ft_freqstatistics(cfg, F{6,:}, F_AB{2,:});
  % (Xpost-ABpost) - (Xpre-ABpre)
  stat{13} = ft_freqstatistics(cfg, Fdiff{7,:}, Fdiff{8,:});
end

save(sprintf(strcat('groupdata', savename, '_%s'),suff), 'stat', 'cmbstat', 'F', 'cmb', 'subj');


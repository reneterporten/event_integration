function rt_corr_mcca(subj, suff)

% Function that correlates prediction values with correlations as output
% fromrt_mcca, resolved over time

datadir = '/project/3012026.13/jansch';
cd(datadir);

predvec = ones(21,1)*-1;
predvec(1) = 1;
predvec(7) = 1;
predvec(12) = 1;
predvec(16) = 1;
predvec(19) = 1;
predvec(21) = 1;

cfg1 = [];
cfg1.appenddim = 'rpt';

cfg2 = [];
cfg2.avgoverrpt = 'no';

d = dir(sprintf('*%s*',[subj suff]));

for k = 1:numel(d)
  
  jobwaste = dir(sprintf('renter*%s*', 'dccn'));
  for waste = 1:numel(jobwaste)
      delete(jobwaste(waste).name)
  end
  
  k
  load(d(k).name);
  fsample = round(1./mean(diff(tlck{1,1}.time))); % assume integer
  n = round(0.25.*fsample);
  for m = 1:size(tlck,1)
    % these are the indivdual pairwise correlations
 
    T{m,k} = ft_selectdata(cfg2, ft_appendtimelock(cfg1, tlck{m,:}));
    if m == 1
        Tpred = T{m,k}.trial;
        Ppred = ones(size(T{m,k}.trial))*predvec(m);
    else
        Tpred = cat(1, Tpred, T{m,k}.trial);
        Ppred = cat(1, Ppred, ones(size(T{m,k}.trial))*predvec(m));
    end
  end 
  
  % Correlating takes a lot of time per channel, maybe use corr functions
  % instead?
  fprintf('Correlating for subject %d\n', k);
  x = reshape(Tpred, 271, 252, 661);
  y = reshape(Ppred, 271, 252, 661);
  for chanx = 1:size(x,1)
    qsubfeval('qsub_trc', x(chanx,:,:), y(chanx,:,:), n, subj, chanx, 'memreq', (1024^3)*24, 'timreq', 60*60);
  end
end


function rt_combine_rsmChan(subj,suff)

datadir = fullfile('/project/3012026.13/jansch/',subj);
cd(datadir);

d = dir(sprintf('sub*%s*',suff));

for k = 1:numel(d)
  k
  load(d(k).name);
  if k == 1
    combChan = zeros(numel(d), numel(corrchan));
  end
  chanstr = extractBetween(d(k).name,"chan_",".mat");
  channum = str2double(chanstr{1});
  combChan(channum,:) = corrchan';
  
end

save(sprintf('combChan_%s',subj), 'combChan');


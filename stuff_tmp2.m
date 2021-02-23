%%
% compute the time-resolved correlations between the paired-repetitions
% the ordering in the 5th dimension is: Apre-Bpre-Xpre - Apst-Bpst-Xpst
clear C;
cliplabel = {'Apre';'Bpre';'Xpre';'Apst';'Bpst';'Xpst'};
for p = 1:nstory
  fprintf('computing correlations for story %d\n', p);
  cnt = 0;
  for k = 1:ncnd
    for m = k:ncnd
      cnt = cnt+1;
      x = reshape(alldat(:,:,p,:,k),nchan, nrpt, ntime);
      y = reshape(alldat(:,:,p,:,m),nchan, nrpt, ntime);
      c = trc(x,y,n);
      
      if k==m
        c(c==1) = nan;
      end
      C(:,:,cnt,p) = squeeze(nanmean(nanmean(c,3),2));
      cmb(cnt,:) = {cliplabel{k} cliplabel{m}};
    end
  end
end

tmp = [];
tmp.label = label;
tmp.time  = data{1}.time;
tmp.dimord = 'chan_time';

tlck = cell(size(C,3), size(C,4));
for k = 1:size(C,3)
  for m = 1:size(C,4)
    tmp.avg = C(:,:,k,m);
    tlck{k,m} = tmp;
  end
end

save_dir = '/project/3012026.13/jansch';
save(fullfile(save_dir, sprintf('%s_mcca_alpha_rsm', subj)), 'tlck', 'cmb');

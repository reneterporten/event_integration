function rt_mcca_FREQ(subj, varargin)

%% 
% get the data
cfg      = ft_getopt(varargin, 'cfg_preproc', []);
root_dir = '/project/3012026.13/processed';
data     = rt_myfreqanalysis(root_dir, subj, 'cfg_preproc', cfg);

%%
% create a 'neighbourhood matrix that determines which channels are to be
% combined
load ctf275_neighb
label = data{1}.label(1:(numel(data{1}.label)/2));
for k = 1:numel(label)
  label{k} = label{k}(1:end-3);
end
[a,b]      = match_str(label, {neighbours.label}');
neighbours = neighbours(b);

pwdir     = pwd;
[~,ftdir] = ft_version;
cd(fullfile(ftdir,'private'));
dummydata.label = label;
cfg = [];
cfg.neighbours = neighbours;
C = channelconnectivity(cfg, dummydata);
C = repmat(C,[1 2]); % channel neighbourhood matrix
cd(pwdir);

%%
% reorganize that data and do mcca
nstory = numel(data)./6;
nfreqs = numel(data{1}.freq);
% ERF data was: trial: [6×542×661 double]
% Freq data is: powspctrm: [6×542×39×45 double], trials, chan, freq, time
% Channels, trials, time*3, nstroy
alldatpre  = zeros(size(C,1), 6, numel(data{1}.time)*3, nfreqs, nstory);
alldatpost = zeros(size(C,1), 6, numel(data{1}.time)*3, nfreqs, nstory);
tic
for k = 1:nstory
  disp(strcat("MCCA Story: ", int2str(k)))
  ix = ((k-1)*3+1):(k*3);
  iy = (ix + nstory*3);
  
  for m = 1:numel(ix)
    %tmpdat = permute(data{ix(m)}.trial, [2 1 3]); % ERF
    tmpdat = permute(data{ix(m)}.powspctrm, [2 1 4 3]); % TFR
    if m==1
      datpre = tmpdat;
    else
      datpre = cat(3, datpre, tmpdat);
    end
  end
  
  for m = 1:numel(iy)
    tmpdat = permute(data{iy(m)}.powspctrm, [2 1 4 3]);
    if m==1
      datpost = tmpdat;
    else
      datpost = cat(3, datpost, tmpdat);
    end
  end
  dat = cat(2, datpre, datpost);
  
  for m = 1:size(C,1)
    tmpdat = dat(C(m,:), :,:,:);
    
    [nchan, nset, ntime, nfreq] = size(tmpdat);
    tmpdat = reshape(tmpdat, nchan*nset, nfreq, ntime);
    for f = 1:nfreq

        tmpdatf = squeeze(tmpdat(:,f,:));
        R = cov(tmpdatf');
        mask = repmat({ones(nchan)}, [1 nset]);
        mask = blkdiag(mask{:});

        S = R.*mask;

        lambda = (trace(S)./size(S,1))./10; % just some value

        % compute the spatial filter and its inverse
        [W, A] = getAW(R, S, 1, ones(1,nset).*nchan, lambda);

        w = mat2cell(squeeze(W), nchan, ones(1,nset));
        w = blkdiag(w{:})';

        tmp = w*tmpdatf;

        alldatpre(m, :, :, f, k)  = tmp(     1:(nset/2), :);
        alldatpost(m, :, :, f, k) = tmp((nset/2+1):nset, :);
    end
  end
end
toc
alldat = cat(3, alldatpre, alldatpost);
clear alldatpre alldatpost

% reorganize the data matrix a bit
[nchan, nrpt, ntime, nfreq, nstory] = size(alldat);
alldat = permute(alldat, [1 2 4 5 3]); % chan, rpt, story, time
alldat = reshape(alldat, [nchan nrpt nfreq nstory ntime/6 6]);

fsample = round(1./mean(diff(data{1}.time))); % assume integer
n = round(0.25.*fsample);

[nchan, nrpt, nfreq, nstory, ntime, ncnd] = size(alldat);

%%
% compute the time-resolved correlations between the paired-repetitions
% the ordering in the 5th dimension is: Apre-Bpre-Xpre - Apst-Bpst-Xpst
clear C;
cliplabel = {'Apre';'Bpre';'Xpre';'Apst';'Bpst';'Xpst'};
for p = 1:nstory
  fprintf('computing correlations for story %d\n', p);
  for k = 1:ncnd

      x = reshape(alldat(:,:,:,p,:,k),nchan, nrpt, nfreq, ntime);

      C(:,:,:,k,p) = squeeze(mean(x, 2));
      cmb(k,:) = {cliplabel{k}};
  end
end

tmp = [];
tmp.label = label;
tmp.time  = data{1}.time;
tmp.dimord = 'chan_freq_time';
tmp.freq = data{1}.freq;

freq = cell(size(C,4), size(C,5));
for k = 1:size(C,4)
  for m = 1:size(C,5)
    tmp.avg = C(:,:,:,k,m);
    freq{k,m} = tmp;
  end
end


save_dir = '/project/3012026.13/jansch';
save(fullfile(save_dir, sprintf('%s_mcca_FREQ', subj)), 'freq', 'cmb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
function [W,A] = getAW(R,S,K,n,lambda)

nset = numel(n);
if isempty(lambda)
  lambda = 0;
end

W = nan+zeros(K,max(n),nset,numel(lambda));
A = nan+zeros(max(n),K,nset,numel(lambda));
R_in = R;
S_in = S;
for m = 1:numel(lambda)
  R = R_in + eye(size(R,1)).*lambda(m);
  S = S_in + eye(size(S,1)).*lambda(m);
  
  % this eigenvalue decomposition gives the unmixing in the columns, so to make
  % it a proper unmixing matrix, to-be-applied to each subject, it should be transposed
  [tempW,~] = eigs((R+R')./2,(S+S')./2,K);
  
  tempW     = normc(tempW);
  tempA     = R*tempW/(tempW'*R*tempW);
  
  sumn = cumsum([0 n(:)']);
  for k = 1:nset
    nchan    = n(k);
    indx     = sumn(k) + (1:nchan);
    W(:,1:nchan,k,m) = (tempW(indx,:))'; %unmixing
    A(1:nchan,:,k,m) = (tempA(indx,:));  %mixing
  end
end

function rho = getrho(R,W,K,n)

for j = 1:size(W,4)
  nset = numel(n);
  
  R(~isfinite(R)) = 0;
  
  tmp = zeros(K*nset,size(R,1));
  sumn = cumsum([0 n(:)']);
  for k = 1:nset
    nchan = n(k);
    for m = 1:K
      tmp((m-1)*nset+k, sumn(k)+(1:nchan)) = W(m,1:nchan,k,j);
    end
  end
  rho(:,:,j) = tmp*R*tmp';
end

function c = trc(x,y,n)

% x = nchan x nrpt x ntime
% y = nchan x nrpt x ntime
%
% computes for each channel the pairwise rpt correlation matrix across
% time, integrating across n time points

[nchanx, nrptx, ntimex] = size(x);
[nchany, nrpty, ntimey] = size(y);
assert(nchanx==nchany);
assert(ntimex==ntimey);

c = nan(nchanx, nrptx, nrpty, ntimex);

x_sum = reshape(ft_preproc_smooth(reshape(x, [], ntimex), n), [nchanx nrptx ntimex]);
y_sum = reshape(ft_preproc_smooth(reshape(y, [], ntimey), n), [nchany nrpty ntimey]);

xsq_sum = reshape(ft_preproc_smooth(reshape(x.^2, [], ntimex), n), [nchanx nrptx ntimex]);
ysq_sum = reshape(ft_preproc_smooth(reshape(y.^2, [], ntimey), n), [nchany nrpty ntimey]);


for k = 1:nrptx
  for m = 1:nrpty
    xy_sum = reshape(ft_preproc_smooth(reshape(x(:,k,:).*y(:,m,:), [], ntimex), n), [nchanx 1 ntimex]);
    
    numer  =  xy_sum(:,1,:)  - (x_sum(:,k,:).*y_sum(:,m,:));
    denomx = (xsq_sum(:,k,:) - (x_sum(:,k,:).^2));
    denomy = (ysq_sum(:,m,:) - (y_sum(:,m,:).^2));

    c(:,k,m,:) = numer./sqrt(denomx.*denomy);
  end
end
 

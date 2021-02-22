function [dataout, tlck] = rt_mcca(data)

% RT_MCCA computes a multiset canonical correlation on an input dataset


%%
% create a 'neighbourhood matrix that determines which channels are to be
% combined for the MCCA computation

% this chunk operates on sensor level data, and assumes the input to be
% planar gradients (uncombined)
% FIXME: for the future it might be relevant to also think of a way to use
% source-level data (parcels).

load ctf275_neighb
label = data.label(1:(numel(data.label)/2));
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
nstory = numel(unique(data.trialinfo(:,4)));
alldatpre  = zeros(size(C,1), 6, numel(data.time{1})*3, nstory);
alldatpost = zeros(size(C,1), 6, numel(data.time{1})*3, nstory);

dataout = keepfields(data, {'trialinfo' 'time'});
for k = 1:numel(data.trial)
  dataout.trial{k} = zeros(size(C,1), numel(data.time{k}));
end

unmixing = zeros([size(C) 12 nstory]); % the 12 is hard coded for now
mixing   = zeros([size(C) 12 nstory]);
for k = 1:nstory
  fprintf('computing weights for story %d/%d\n', k, nstory);
  
  % identify the epochs that belong, for a given story to clips A/B/X for
  % the pre and post linking event 
  ix = find(data.trialinfo(:,4)==k & ismember(data.trialinfo(:,2), [1 2 3]));
  iy = find(data.trialinfo(:,4)==k & ismember(data.trialinfo(:,2), [5 6 7]));
  
  % hard coded for now, works with the current data organization
  ix = reshape(ix, [], 3);
  iy = reshape(iy, [], 3);
  
  [nchan, ntime] = size(data.trial{1});
  datpre = zeros(nchan, size(ix,1), ntime*size(ix,2));
  for m = 1:size(ix,1)
    for mm = 1:size(ix,2)
      datpre(:,m,(mm-1)*ntime+(1:ntime)) = data.trial{ix(m,mm)} - nanmean(data.trial{ix(m,mm)},2);
    end
  end
  
  [nchan, ntime] = size(data.trial{1});
  datpost = zeros(nchan, size(iy,1), ntime*size(iy,2));
  for m = 1:size(iy,1)
    for mm = 1:size(iy,2)
      datpost(:,m,(mm-1)*ntime+(1:ntime)) = data.trial{iy(m,mm)} - nanmean(data.trial{iy(m,mm)},2);
    end
  end
  dat = cat(2, datpre, datpost);
  
  for m = 1:size(C,1)
    tmpdat = dat(C(m,:), :,:);
    
    [nchan, nset, ntime] = size(tmpdat);
    tmpdat = reshape(tmpdat, nchan*nset, ntime);
    
    R = cov(tmpdat');
    mask = repmat({ones(nchan)}, [1 nset]);
    mask = blkdiag(mask{:});
    
    S = R.*mask;

    lambda = (trace(S)./size(S,1))./10; % just some value
    
    % compute the spatial filter and its inverse
    [W, A] = getAW(R, S, 1, ones(1,nset).*nchan, lambda);
    
    unmixing(m,C(m,:),:,k) = squeeze(W);
    mixing(m,C(m,:),:,k)   = squeeze(A);
    
    %w = mat2cell(squeeze(W), nchan, ones(1,nset));
    %w = blkdiag(w{:})';
    
    %a = mat2cell(squeeze(A), nchan, ones(1,nset));
    %a = blkdiag(a{:});
  end
end

% here, try and align the polarity of the weights such that 
% 1) across stories the signals are polarity consistent, and
% 2) across neighbouring channels the signals are polarity consistent
%
% first try to align the polarity per story across channels
mixing_flipped = mixing;
unmixing_flipped = unmixing;
for k = 1:nstory
  fprintf('aligning polarity across channels for story %d/%d\n', k, nstory);
  input = reshape(mixing(:,:,:,k), size(mixing,1), []);
  [output, flipped] = polarity_align(input, 0); 
  
  siz = size(mixing);
  mixing_flipped(:,:,:,k) = reshape(output, [size(output,1), siz(2:3)]);
  unmixing_flipped(:,:,:,k) = repmat(flipped, [1 siz(2:3)]).*unmixing(:,:,:,k);
end

% combine the weights with the corresponding input data.
index = cell(1,nstory);
for k = 1:nstory
  fprintf('computing canonical components for story %d/%d\n', k, nstory);
  
  % identify the epochs that belong, for a given story to clips A/B/X for
  % the pre and post linking event 
  ix = find(data.trialinfo(:,4)==k & ismember(data.trialinfo(:,2), [1 2 3]));
  iy = find(data.trialinfo(:,4)==k & ismember(data.trialinfo(:,2), [5 6 7]));
  
  % hard coded for now, works with the current data organization
  ix = reshape(ix, [], 3);
  iy = reshape(iy, [], 3);
  
  for m = 1:size(ix,1)
    for mm = 1:size(ix,2)
      dataout.trial{ix(m,mm)} = unmixing_flipped(:,:,m,mm) * data.trial{ix(m,mm)};
      dataout.trial{iy(m,mm)} = unmixing_flipped(:,:,m,mm) * data.trial{iy(m,mm)};
     
      % facilitates covariance estimate for final polarity check/adjustment
      dataout.trial{ix(m,mm)} = dataout.trial{ix(m,mm)} - nanmean(dataout.trial{ix(m,mm)},2);
      dataout.trial{iy(m,mm)} = dataout.trial{iy(m,mm)} - nanmean(dataout.trial{iy(m,mm)},2);
    end
  end
  
  index{k} = [ix(:);iy(:)];
end

% compute the components average per story, and check whether the
% corresponding signals are aligned
% dat = zeros([size(dataout.trial{1}) nstory]);
% for k = 1:nstory
%   dat(:,:,k) = mean(cat(3, dataout.trial{index{k}}),3);
% end
% 
% for k = 1:nstory
%   for m = 1:nstory
%     if k~=m
%     cdat(:,k,m) = sum(dat(:,:,k).*dat(:,:,m),2)./sqrt(sum(dat(:,:,k).*dat(:,:,k),2) .* sum(dat(:,:,m).*dat(:,:,m),2));
%     end
%   end
% end

dataout.label = label;

% compute the average across clips, that's basically for free
conds = [1 2 3 5 6 7];
for k = 1:numel(conds)
  tmpcfg        = [];
  tmpcfg.trials = find(dataout.trialinfo(:,2)==conds(k));
  tmpcfg.preproc.baselinewindow = [-0.1 0];
  tmpcfg.preproc.demean         = 'yes';
  tlck(k) = ft_timelockanalysis(tmpcfg, dataout);
end

for k = 1:6
  tlck(k).time = tlck(1).time;
end

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
 
function [output, flipped] = polarity_align(input, doplot)

% input is a nchan*something matrix, and the idea is to adjust the polarity
% of individual channels such, that the entries in the covariance matrix
% that are non-zero are as aligned as possible, as reflected in positive
% values

if nargin<2
  doplot = false;
end

ok = false;
output = input;
sprevious = zeros(1,size(input,1));

sC = get_sC(input);

err = zeros(0,1);
flipped = ones(size(input,1),1);

while numel(err)<5 || sum(err(end-4:end))>2 
  [m, sel] = min(sC);
  output(sel,:) = -output(sel,:);
  flipped(sel)  = -flipped(sel);
  
  sC = get_sC(output);
  
  if doplot
    hold off;
    plot(sC)
    hold on;
    plot(sprevious);drawnow
  end
  
  %pause;
  err(end+1) = sum(sC<0);
  
  sprevious = sC;
  
end

function sC = get_sC(input)

C = input*input';
C = C-diag(diag(C));
sC = sum(sign(C))./sum(sign(C)~=0);

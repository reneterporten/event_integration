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
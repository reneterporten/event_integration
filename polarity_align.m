function [output, flipped] = polarity_align(input, doplot)

% input is a nchan*something matrix, and the idea is to adjust the polarity
% of individual channels such, that the entries in the covariance matrix
% that are non-zero are as aligned as possible, as reflected in positive
% values

if nargin<1
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

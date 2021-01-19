function [subtrdata] = rt_subtrmean(data)

% Function to subtract the channel specific mean from the data
% Potentially accounts from differences in amplitude offset as a
% consequence of low frequency oscillations

data = reshape(data, size(data,2), size(data,1), size(data, 3));
subtrdata = zeros(size(data));
for chan = 1:size(data,1)
    
    chandata            = squeeze(data(chan, :, :));
    chandata            = chandata - repmat(mean(chandata,2),1,size(chandata,2));
    subtrdata(chan,:,:) = chandata;
    
end
clear data
subtrdata = reshape(subtrdata, size(subtrdata,2), size(subtrdata,1), size(subtrdata, 3));
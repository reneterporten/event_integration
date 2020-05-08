function [predictionRDMxlsx] = rt_predictionRDMxlsx(predata, postdata)

% Function to read in a template prediction RDM from excel file
% Applies the tempalte RDM to overall data structure that fits the neural
% data

% Function parameters are customized to the study

predictionRDMxlsx = NaN(size(predata,1)*24, size(predata,1)*24);

for preRDMs = 1:18:(size(predictionRDMxlsx,1)/2)
    
    predictionRDMxlsx(preRDMs:(preRDMs+17), preRDMs:(preRDMs+17)) = predata;
    
end

for postRDMs = 217:18:(size(predictionRDMxlsx,1))
    
    predictionRDMxlsx(postRDMs:(postRDMs+17), postRDMs:(postRDMs+17)) = postdata;
    
end
function [predictionRDMxlsx] = rt_predictionRDMsource(cfg, predata, postdata)

% Function to read in a template prediction RDM from excel file
% Applies the tempalte RDM to overall data structure that fits the neural
% data

% Function parameters are customized to the study

predictionRDMxlsx = NaN(size(predata,1)*24, size(predata,1)*24);

for preRDMs = 1:3:(size(predictionRDMxlsx,1)/2)
    
    predictionRDMxlsx(preRDMs:(preRDMs+2), preRDMs:(preRDMs+2)) = predata;
    
end

for postRDMs = 37:3:(size(predictionRDMxlsx,1))
    
    predictionRDMxlsx(postRDMs:(postRDMs+2), postRDMs:(postRDMs+2)) = postdata;
    
end

if cfg.offdiag == true
    % Apply off diagional (lower) - diagonal
    offDiagIdx = 1:3:(size(predictionRDMxlsx,1)/2);
    offDiagCount = 1;
    for offdiag = 37:3:(size(predictionRDMxlsx,1))

        predictionRDMxlsx(offdiag:(offdiag+2), offDiagIdx(offDiagCount):(offDiagIdx(offDiagCount)+2)) = postdata;
        offDiagCount = offDiagCount + 1;

    end

    % Apply off diagional (upper) - diagonal
    offDiagIdx = 1:3:(size(predictionRDMxlsx,1)/2);
    offDiagCount = 1;
    for offdiag = 37:2:(size(predictionRDMxlsx,1))

        predictionRDMxlsx(offDiagIdx(offDiagCount):(offDiagIdx(offDiagCount)+17), offdiag:(offdiag+17)) = postdata;
        offDiagCount = offDiagCount + 1;

    end
end
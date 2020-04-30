function [predictionRDM] = rt_predictionRDM(pred_template, all_data)

% Function that returns a prediction RDM structure that matches the trial
% structure per spatial & temporal searchlight in the data

% all_data should be in format:
% all_data = {erp_A_pre; erp_B_pre; erp_C_pre};

% pred_template should be in format:
% pred_template = [NaN, -1, 1; -1, NaN, 1; 1, 1, NaN];

allTrialNr      = {};
for dataNr = 1:size(all_data, 1)
    allTrialNr{dataNr} = size(all_data{dataNr}.trialinfo,1);
end

trial_vec_col   = [];
for trialvecNr = 1:size(allTrialNr, 2)
    trial_vec_col = [trial_vec_col, allTrialNr{trialvecNr}];
end

trial_vec_row   = [];
trial_vec_row   = trial_vec_col';

predictionStruc = {};
predictionStruc = cell(size(pred_template));

for rowStruc = 1:size(predictionStruc,1)
    for colStruc = 1:size(predictionStruc,2)
        predictionStruc{rowStruc, colStruc} = [trial_vec_row(rowStruc, 1), trial_vec_col(1, colStruc)];
    end
end

for rowStruc2 = 1:size(predictionStruc,1)
    for colStruc2 = 1:size(predictionStruc,2)
        predictionStruc{rowStruc2, colStruc2} = ones(predictionStruc{rowStruc2, colStruc2}).*pred_template(rowStruc2, colStruc2);
    end
end

predictionRDM       = [];
predictionRDM_row   = [];
for rowStruc3 = 1:size(predictionStruc,1)
    for colStruc3 = 1:size(predictionStruc,2)
        predictionRDM_row = [predictionRDM_row, predictionStruc{rowStruc3, colStruc3}];
    end
    predictionRDM       = [predictionRDM; predictionRDM_row];
    predictionRDM_row   = [];
end

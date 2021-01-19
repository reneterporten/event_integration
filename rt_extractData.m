function extractedData = rt_extractData(exdata)

idxDiag         = logical(eye(size(exdata)));
exdata(idxDiag) = NaN; % Replacing man diagonal elements with NaN

% Select the pre, post and prepost data from initial structure
datapoints  = length(exdata);
predata     = exdata(1:datapoints/2, 1:datapoints/2);
postdata    = exdata((datapoints/2)+1:datapoints, (datapoints/2)+1:datapoints);
prepostdata = exdata((datapoints/2)+1:datapoints, 1:datapoints/2);

% Create identifier to target comparisons of conditions (fo each phase condition)
dummyMatrix = [ones(6,6)*11, ones(6,6)*12, ones(6,6)*13; 
    ones(6,6)*21, ones(6,6)*22, ones(6,6)*23; 
    ones(6,6)*31, ones(6,6)*32, ones(6,6)*33];
dummyMatrixSmall = [11, 12, 13; 
    21, 22, 23; 
    31, 32, 33];
idxMatrix = repmat(dummyMatrix,[12,12]);

% Adjust idxMatrix to regard only within story comparisons
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'narrativeInsight.mat'), 'narrativeInsight')
narrativeInsight = narrativeInsight(1:datapoints/2, 1:datapoints/2);
celcounter = 1;
for narRow = 1:(length(narrativeInsight)/6)
    narrativeInsight(celcounter:(celcounter+5), celcounter:(celcounter+5)) = 1;
    celcounter = celcounter + 6;
end
idxMatrix(isnan(narrativeInsight)) = NaN;

% Condition pre
% Select only lower triangle of data
allIdx = cell(6,1);
allIdx{1} = ismember(idxMatrix, [11]); % AA
allIdx{2} = ismember(idxMatrix, [21]); % AB
allIdx{3} = ismember(idxMatrix, [31]); % AX
allIdx{4} = ismember(idxMatrix, [33]); % XX
allIdx{5} = ismember(idxMatrix, [22]); % BB
allIdx{6} = ismember(idxMatrix, [32]); % BX

% Prepare structures that contain the mean of the cell comparisons
prestruc        = ones(6,1);
poststruc       = ones(6,1);
prepoststruc    = ones(6,1);

% Loop through structure and calculate mean
for struc = 1:length(prestruc)
    prestruc(struc,1)       = nanmean(predata(allIdx{struc}));
    poststruc(struc,1)      = nanmean(postdata(allIdx{struc}));
    prepoststruc(struc,1)   = nanmean(prepostdata(allIdx{struc}));
end

pre_extr_mat        = zeros(3,3);
post_extr_mat       = zeros(3,3);
prepost_extr_mat    = zeros(3,3);
comborder = [11, 21, 31, 33, 22, 32]';
for sortdata = 1:length(comborder)
    currentNumber = comborder(sortdata);
    pre_extr_mat(ismember(dummyMatrixSmall, currentNumber))     = prestruc(sortdata,1);
    post_extr_mat(ismember(dummyMatrixSmall, currentNumber))    = poststruc(sortdata,1);
    prepost_extr_mat(ismember(dummyMatrixSmall, currentNumber)) = prepoststruc(sortdata,1);
end

% Mirror lower triangle values to upper triangle
pre_extr_mat        = (pre_extr_mat+pre_extr_mat') - eye(size(pre_extr_mat,1)).*diag(pre_extr_mat);
post_extr_mat       = (post_extr_mat+post_extr_mat') - eye(size(post_extr_mat,1)).*diag(post_extr_mat);
prepost_extr_mat    = (prepost_extr_mat+prepost_extr_mat') - eye(size(prepost_extr_mat,1)).*diag(prepost_extr_mat);

% Data structure resulting from extraction
extractedData = [pre_extr_mat, prepost_extr_mat; prepost_extr_mat, post_extr_mat];






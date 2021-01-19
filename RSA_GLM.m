%% Using extracted data for GLM analyses
%% Default path

addpath /home/common/matlab/fieldtrip/
addpath /home/common/matlab/fieldtrip/qsub/
addpath /project/3012026.13/scripts_RT/
ft_defaults

root_dir     = '/project/3012026.13/processed/';
save_dir     = '/project/3012026.13/processed_RT/time_lock/';

subjects = {'sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-016',...
    'sub-017','sub-018','sub-019','sub-020','sub-022','sub-023',...
    'sub-024','sub-025','sub-026','sub-027','sub-028','sub-029','sub-031'...
    'sub-032','sub-033','sub-034','sub-035','sub-036','sub-037','sub-038'};


%% Load in neural data

% This helps to select channels
load(fullfile('/project/3012026.13/processed_RT/', 'template_timelock.mat'));
cfg                 = [];
cfg.trials          = 'all';
cfg.avgoverrpt      = 'yes';
cfg.avgovertime     = 'no';
template_data       = ft_selectdata(cfg, template_timelock);
template_data.avg = template_data.trial;
clear template_timelock

% Create data structure that contains RSA data of all subjects
% Allocate overall structure
ignSub = [4, 10, 12, 20, 30, 32];
allDataGLM = cell((length(subjects)-length(ignSub)),1);

rsa_count = 1;
for subData = 1:length(subjects)
    
    if ~ismember(subData, ignSub)

        disp(strcat('Creating GLM structure:', int2str(subData)))

        subj = subjects{subData};

        load(fullfile('/project/3012026.13/processed_RT/RSA Results/', subj,'resulting_data_GLM.mat'));

        allDataGLM{rsa_count}   = resulting_data_GLM.searchlight;
        rsa_count               = rsa_count + 1;
    
    end

end

neuralSim_avg = zeros(size(allDataGLM{1}));
for GLMchan = 1:size(allDataGLM{1},1)
    for GLMtime = 1:size(allDataGLM{1},2)
        for GLMsubs = 1:length(allDataGLM)
            neuralSim_avg(GLMchan, GLMtime,:,:) = neuralSim_avg(GLMchan, GLMtime,:,:) + allDataGLM{GLMsubs}(GLMchan, GLMtime,:,:);
        end
        neuralSim_avg(GLMchan, GLMtime,:,:) = neuralSim_avg(GLMchan, GLMtime,:,:)./length(allDataGLM);
    end
end

%% Prepare structure for extraction for R

% Select channels
antchan     = {'MLC12', 'MLC13', 'MLC14', 'MLC15', 'MLC21', 'MLC22', 'MLC23', 'MLC24', 'MLC31', 'MLC41', 'MLC51', 'MLC52', 'MLC53', 'MLC61', 'MLC62', 'MLF23', 'MLF31', 'MLF32', 'MLF41', 'MLF42', 'MLF43', 'MLF51', 'MLF52', 'MLF53', 'MLF61', 'MLF62', 'MLF63', 'MRC11', 'MRC12', 'MRC13', 'MRC21', 'MRC22', 'MRC23', 'MRC31', 'MRC41', 'MRC51', 'MRC52', 'MRC53', 'MRC61', 'MRC62', 'MRF31', 'MRF32', 'MRF41', 'MRF42', 'MRF43', 'MRF51', 'MRF52', 'MRF61', 'MRF62', 'MZC01', 'MZC02', 'MZC03', 'MZF02', 'MZF03'};
mychannels  = antchan;
channels    = cell2mat(template_data.label);

data4table = [];
for Rrunner = 1:length(allDataGLM)
    Rdataselect = allDataGLM{Rrunner}(ismember(channels, mychannels),:,:,:); % Select data from channels
    Rdataselect = squeeze(mean(Rdataselect,1)); % Calculate mean over selected channels
    Rdataselect = squeeze(mean(Rdataselect,1)); % Calculate mean of over the 15 time windows

    AB_pre  = Rdataselect(1,2);
    AX_pre  = Rdataselect(1,3);
    AB_post = Rdataselect(4,5);
    AX_post = Rdataselect(4,6);
    subject = Rrunner;
    data4table = [data4table; Rrunner, AB_pre, AX_pre, AB_post, AX_post];
end

% Create table that can be interpreted by R
data4GLM = array2table(data4table,'VariableNames',{'Subject','AB_pre','AX_pre','AB_post','AX_post'});
% Save table
writetable(data4GLM,fullfile('/project/3012026.13/processed_RT/RSA Results/Data for GLM/data4GLM.txt'))

%%

% Devide the brain into frontal,L/R temporal and posterior channels
% Plot RSMs per time window

antchan = {'MLC12', 'MLC13', 'MLC14', 'MLC15', 'MLC21', 'MLC22', 'MLC23', 'MLC24', 'MLC31', 'MLC41', 'MLC51', 'MLC52', 'MLC53', 'MLC61', 'MLC62', 'MLF23', 'MLF31', 'MLF32', 'MLF41', 'MLF42', 'MLF43', 'MLF51', 'MLF52', 'MLF53', 'MLF61', 'MLF62', 'MLF63', 'MRC11', 'MRC12', 'MRC13', 'MRC21', 'MRC22', 'MRC23', 'MRC31', 'MRC41', 'MRC51', 'MRC52', 'MRC53', 'MRC61', 'MRC62', 'MRF31', 'MRF32', 'MRF41', 'MRF42', 'MRF43', 'MRF51', 'MRF52', 'MRF61', 'MRF62', 'MZC01', 'MZC02', 'MZC03', 'MZF02', 'MZF03'};
postchan = {'MLO11', 'MLO12', 'MLO13', 'MLO21', 'MLO22', 'MLO23', 'MLO24', 'MLO31', 'MLO32', 'MLO41', 'MLO42', 'MLO43', 'MLP21', 'MLP31', 'MLP32', 'MLP33', 'MLP41', 'MLP42', 'MLP43', 'MLP51', 'MLP52', 'MLP53', 'MLP54', 'MRO11', 'MRO12', 'MRO13', 'MRO14', 'MRO21', 'MRO22', 'MRO23', 'MRO24', 'MRO31', 'MRO32', 'MRO33', 'MRO41', 'MRO42', 'MRO43', 'MRO44', 'MRO53', 'MRP21', 'MRP31', 'MRP32', 'MRP33', 'MRP41', 'MRP42', 'MRP43', 'MRP44', 'MRP51', 'MRP52', 'MRP53', 'MRP54', 'MRP55', 'MZO01', 'MZO02', 'MZP01'};
leftchan = {'MLC15', 'MLC16', 'MLC17', 'MLC24', 'MLC25', 'MLF35', 'MLF44', 'MLF45', 'MLF46', 'MLF54', 'MLF55', 'MLF56', 'MLF63', 'MLF64', 'MLF65', 'MLF66', 'MLF67', 'MLO14', 'MLP35', 'MLP43', 'MLP44', 'MLP45', 'MLP54', 'MLP55', 'MLP56', 'MLP57', 'MLT11', 'MLT12', 'MLT13', 'MLT14', 'MLT15', 'MLT16', 'MLT21', 'MLT22', 'MLT23', 'MLT24', 'MLT25', 'MLT26', 'MLT31', 'MLT32', 'MLT33', 'MLT34', 'MLT35', 'MLT36', 'MLT41', 'MLT42', 'MLT43', 'MLT44', 'MLT45', 'MLT53'};
rightchan = {'MRC13', 'MRC14', 'MRC15', 'MRC16', 'MRC17', 'MRC22', 'MRC23', 'MRC24', 'MRC25', 'MRC31', 'MRC32', 'MRC42', 'MRF25', 'MRF34', 'MRF35', 'MRF44', 'MRF45', 'MRF46', 'MRF52', 'MRF53', 'MRF54', 'MRF55', 'MRF56', 'MRF62', 'MRF63', 'MRF64', 'MRF65', 'MRF67', 'MRO14', 'MRP23', 'MRP33', 'MRP34', 'MRP35', 'MRP42', 'MRP43', 'MRP44', 'MRP45', 'MRP53', 'MRP54', 'MRP55', 'MRP56', 'MRP57', 'MRT11', 'MRT12', 'MRT13', 'MRT14', 'MRT15', 'MRT16', 'MRT21', 'MRT22', 'MRT23', 'MRT24', 'MRT25', 'MRT26', 'MRT27', 'MRT31', 'MRT32', 'MRT33', 'MRT34', 'MRT35', 'MRT36', 'MRT42', 'MRT43', 'MRT44', 'MRT45', 'MRT53', 'MRT54'};

% Anterior channels
mychannels  = antchan;
channels    = cell2mat(template_data.label);
dataselect  = neuralSim_avg(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    caxis([-0.25 0.25]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('Anterior Electrodes Neural RSMs')

% Posterior channels
mychannels  = postchan;
channels    = cell2mat(template_data.label);
dataselect  = neuralSim_avg(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    caxis([-0.25 0.25]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('Posterior Electrodes Neural RSMs')

% left channels
mychannels  = leftchan;
channels    = cell2mat(template_data.label);
dataselect  = neuralSim_avg(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    caxis([-0.25 0.25]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('Left Lateral Electrodes Neural RSMs')

% left channels
mychannels  = rightchan;
channels    = cell2mat(template_data.label);
dataselect  = neuralSim_avg(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    caxis([-0.25 0.25]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('Right Lateral Electrodes Neural RSMs')


%% Plot the difference between cells (pre-post)

% Calculate difference

neuralSim_post_pre_diff = neuralSim_avg(:,:,4:6,4:6) - neuralSim_avg(:,:,1:3,1:3);

% Anterior channels
mychannels  = antchan;
channels    = cell2mat(template_data.label);
dataselect  = neuralSim_post_pre_diff(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

% Barplot
figure
for images = 1:size(dataselect,1)
    ax(images)  = subplot(5,3,images);
    withinEvent = diag(squeeze(dataselect(images,:,:)))'; % AA, BB, XX
    acrossEvent = vectorizeSimmat(squeeze(dataselect(images,:,:))); % AB, AX, BX
    allComp     = [withinEvent, acrossEvent];
    bar(allComp)
    set(gca, 'XTickLabel', {'AA' 'BB' 'XX' 'AB' 'AX' 'BX'})
    set(gcf, 'Position',  [100, 100, 500, 800])
    xtickangle(45)
    ylim([-9*10^(-3) 9*10^(-3)])
    title(strcat('T', num2str(images)))
end
sgtitle('Anterior Electrodes Post-Pre Phase Similarity')

% Anterior channels
mychannels  = postchan;
channels    = cell2mat(template_data.label);
dataselect  = neuralSim_post_pre_diff(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

% Barplot
figure
for images = 1:size(dataselect,1)
    ax(images)  = subplot(5,3,images);
    withinEvent = diag(squeeze(dataselect(images,:,:)))'; % AA, BB, XX
    acrossEvent = vectorizeSimmat(squeeze(dataselect(images,:,:))); % AB, AX, BX
    allComp     = [withinEvent, acrossEvent];
    bar(allComp)
    set(gca, 'XTickLabel', {'AA' 'BB' 'XX' 'AB' 'AX' 'BX'})
    set(gcf, 'Position',  [100, 100, 500, 800])
    xtickangle(45)
    ylim([-9*10^(-3) 9*10^(-3)])
    title(strcat('T', num2str(images)))
end
sgtitle('Posterior Electrodes Post-Pre Phase Similarity')

% Anterior channels
mychannels  = leftchan;
channels    = cell2mat(template_data.label);
dataselect  = neuralSim_post_pre_diff(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

% Barplot
figure
for images = 1:size(dataselect,1)
    ax(images)  = subplot(5,3,images);
    withinEvent = diag(squeeze(dataselect(images,:,:)))'; % AA, BB, XX
    acrossEvent = vectorizeSimmat(squeeze(dataselect(images,:,:))); % AB, AX, BX
    allComp     = [withinEvent, acrossEvent];
    bar(allComp)
    set(gca, 'XTickLabel', {'AA' 'BB' 'XX' 'AB' 'AX' 'BX'})
    set(gcf, 'Position',  [100, 100, 500, 800])
    xtickangle(45)
    ylim([-9*10^(-3) 9*10^(-3)])
    title(strcat('T', num2str(images)))
end
sgtitle('Left Lateral Electrodes Post-Pre Phase Similarity')

% Anterior channels
mychannels  = rightchan;
channels    = cell2mat(template_data.label);
dataselect  = neuralSim_post_pre_diff(ismember(channels, mychannels),:,:,:);
dataselect  = squeeze(mean(dataselect, 1));

% Barplot
figure
for images = 1:size(dataselect,1)
    ax(images)  = subplot(5,3,images);
    withinEvent = diag(squeeze(dataselect(images,:,:)))'; % AA, BB, XX
    acrossEvent = vectorizeSimmat(squeeze(dataselect(images,:,:))); % AB, AX, BX
    allComp     = [withinEvent, acrossEvent];
    bar(allComp)
    set(gca, 'XTickLabel', {'AA' 'BB' 'XX' 'AB' 'AX' 'BX'})
    set(gcf, 'Position',  [100, 100, 500, 800])
    xtickangle(45)
    ylim([-9*10^(-3) 9*10^(-3)])
    title(strcat('T', num2str(images)))
end
sgtitle('Right Lateral Electrodes Post-Pre Phase Similarity')


%%
% Imagesc
figure
for images = 1:size(dataselect,1)
    ax(images) = subplot(5,3,images);
    imagesc(squeeze(dataselect(images,:,:)))
    %caxis([-0.02 0.02]);
    set(gcf, 'Position',  [100, 100, 500, 800])
    title(strcat('T', num2str(images)))
end
sgtitle('Anterior Electrodes Post-Pre Phase Similarity')


%%
% Check some channels
cfg                 = [];
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'avg';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
ft_multiplotER(cfg, template_data)

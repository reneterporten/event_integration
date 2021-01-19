%% Pipeline to create all prediction RSMs
%% Load default data

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


%% Visual similarity 1 (within vs. across event comparisons)

sanitydata      = xlsread('/project/3012026.13/scripts_RT/sanity_pre.xlsx'); % Visual similarity
cfg             = [];
cfg.offdiag     = true; % Use postdata structure for off diagional (phase comparison)
visual_similarity_1   = rt_predictionRDMxlsx(cfg, sanitydata, sanitydata);
visual_similarity_1(isnan(visual_similarity_1))             = -1; % Use -1 instead of NaN for visual similarity
visual_similarity_1(find(eye(size(visual_similarity_1))))   = NaN;% Main diagonal same item comparisons are NaN

imagesc(visual_similarity_1)

save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'visual_similarity_1.mat'), 'visual_similarity_1')
clear visual_similarity_1

% Plot
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'visual_similarity_1.mat'))
[nr,nc] = size(visual_similarity_1);

pcolor([visual_similarity_1 nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
title('Within vs across event similarity')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 


%% Visiual similarity 2 (within vs. across event comparisons, per subject video data)

errorVisSubs = [];
for vis2sub = 1:length(subjects)
    
    try
        disp(strcat('Visual Similarity 2: ', subjects{vis2sub}))

        subj = subjects{vis2sub};
        load(fullfile(root_dir, [subj,'_dataclean.mat']));

        [vs2_magn, vs2_phas] = rt_predictionRDM_videodata(dataclean);

        save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_vs2_magn.mat']), 'vs2_magn')
        save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_vs2_phas.mat']), 'vs2_phas')
        clear vs2_magn
        clear vs2_phas
        clear dataclean    
    catch        
        disp(strcat('Error: ', subjects{vis2sub}))
        errorVisSubs = [errorVisSubs, vis2sub];
        % It appears that following vis2sub have inconsistent trial
        % numbers: 4    10    12    20    30    32
    end

end

% Plot magn
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', ['sub-004', '_vs2_magn.mat']))
[nr,nc] = size(vs2_magn);

pcolor([vs2_magn nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
title('Magnitude')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu')))

% Plot phase
% Plot
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', ['sub-004', '_vs2_phas.mat']))
[nr,nc] = size(vs2_phas);

pcolor([vs2_phas nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
title('Phase')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu'))) 


%% Narrative insight 1 (interaction)

predataNoBX     = xlsread('/project/3012026.13/scripts_RT/prediction_pre_noBX.xlsx');
postdataNoBX    = xlsread('/project/3012026.13/scripts_RT/prediction_post_noBX.xlsx');
predata         = xlsread('/project/3012026.13/scripts_RT/prediction_pre.xlsx');
postdata        = xlsread('/project/3012026.13/scripts_RT/prediction_post.xlsx');
%Create prediction RDM that matches trial x trial structure of data
cfg             = [];
cfg.offdiag     = false; % Use postdata structure for off diagional (phase comparison)
narrativeInsight        = rt_predictionRDMxlsx(cfg, predata, postdata);
narrativeInsight_NoBX   = rt_predictionRDMxlsx(cfg, predataNoBX, postdataNoBX);

save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'narrativeInsight.mat'), 'narrativeInsight')
save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'narrativeInsight_NoBX.mat'), 'narrativeInsight_NoBX')

% Plot
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'narrativeInsight.mat'))
[nr,nc] = size(narrativeInsight);

pcolor([narrativeInsight nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
title('NarrativeInsight Interaction')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu')))


%% Narrative insight 2 (behavior per subject)

errorBehav = [];
for behavsub = 1:length(subjects)
    
    try
        disp(strcat('Behavior insight: ', subjects{behavsub}))

        subj            = subjects{behavsub};
        subj            = strcat(subj(1:3), subj(5:length(subj)));
        sub_behav_path  = fullfile('/project/3012026.13/processed_RT/logfiles/', [subj,'_Relatedness_answers.txt']);
        behaviorInsight = rt_predictionRDMbehav(sub_behav_path);

        % Prepare No BX structure
        nanStruc                = nan(18,18);
        nanStruc(7:12,13:18)    = 1;
        nanStruc(13:18,7:12)    = 1;
        baseStruc               = ones(18,18);
        baseStruc(nanStruc == 1) = NaN;

        % Apply no BX structure
        behaviorInsight_NoBX = nan(length(behaviorInsight), length(behaviorInsight));
        for baseRow = 1:18:length(behaviorInsight)
            disp(baseRow)
            for baseCol = 1:18:length(behaviorInsight)
                behaviorInsight_NoBX(baseRow:baseRow+17, baseCol:baseCol+17) = behaviorInsight(baseRow:baseRow+17, baseCol:baseCol+17).*baseStruc;
            end
        end

        subj            = subjects{behavsub};
        save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_behaviorInsight.mat']), 'behaviorInsight')
        save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_behaviorInsight_NoBX.mat']), 'behaviorInsight_NoBX')
    catch
        disp(strcat('Error: ', subjects{behavsub}))
        errorBehav = [errorBehav, behavsub];
    end

end

% Plot
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', ['sub-004', '_behaviorInsight.mat']))
[nr,nc] = size(behaviorInsight);

pcolor([behaviorInsight nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
title('Narrative Insight Behavior')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu')))


%% Narrative mismatch

errorMisSub = [];
for misSub = 1:length(subjects)
    
    try
        disp(strcat('Narrative mismatch: ', subjects{misSub}))

        subj = subjects{misSub};
        load(fullfile(root_dir, [subj,'_dataclean.mat']));

        narrativeMismatch = rt_predictionRDM_mismatch(dataclean);

        save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', [subj, '_narrativeMismatch.mat']), 'narrativeMismatch')
        clear narrativeMismatch
        clear dataclean
    catch
        disp(strcat('Error: ', subjects{misSub}))
        errorMisSub = [errorMisSub, misSub];
    end

end

% Plot
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', ['sub-004', '_narrativeMismatch.mat']))
[nr,nc] = size(narrativeMismatch);

pcolor([narrativeMismatch nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
title('Narrative Mismatch')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu')))


%% (Un)linking mechanism

un_linkMech = rt_predictionRDM_linkmech();
save(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'un_linkMech.mat'), 'un_linkMech')

% Plot
load(fullfile('/project/3012026.13/scripts_RT/Prediction RSMs', 'un_linkMech.mat'))
[nr,nc] = size(un_linkMech);

pcolor([un_linkMech nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
title('(Un)linking Mechanism')
colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);
colormap(flipud(brewermap(64,'RdBu')))

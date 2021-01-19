%% Analyse ERPs to investigate surprisal
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

%%
tic
mytimedata  = load(fullfile(root_dir, ['sub-004','_dataclean.mat']));
mytime      = mytimedata.dataclean.time{1};
clear mytimedata 

allTrialsPre = [];
allTrialsPost = [];
for iSubject = 1:length(subjects)
    
    disp(strcat('Calculating PP:', int2str(iSubject)))
    
    % Load Data
    subj                = subjects{iSubject};
    
    % load data, remove bad trials and demean & detrend
    load(fullfile(root_dir, [subj,'_dataclean.mat']));
    
    cfg             = [];
    cfg.hpfilter    = 'yes';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 35;
    cfg.hpfreq      = 0.2;
    dataclean       = ft_preprocessing(cfg, dataclean);
    
    mytimestruc = cell(1,length(dataclean.trialinfo));
    for timings = 1:length(mytimestruc)
        mytimestruc{timings} = mytime;
    end
    
    cfg             = [];
    cfg.time        = mytimestruc;
    dataclean       = ft_resampledata(cfg, dataclean);

    % Recode trialinfo for linking events
    indexLink   = 1:18:length(dataclean.trialinfo);
    do_dont = -1;
    for indices = 1:24  %Number of stories * 2
        if do_dont == 1
            dataclean.trialinfo(indexLink(indices):indexLink(indices)+5, 2) = 4;
            dataclean.trialinfo(indexLink(indices):indexLink(indices)+5, 3) = 4;
            dataclean.trialinfo(indexLink(indices):indexLink(indices)+5, 5) = 3;
            indexLink = indexLink+6;
        end
        do_dont = do_dont*-1;
    end

    cfg         = [];
    cfg.trials  = find(dataclean.trialinfo(:,5)==1);
    data_pre    = ft_selectdata(cfg, dataclean);

    cfg         = [];
    cfg.trials  = find(dataclean.trialinfo(:,5)==2);
    data_post   = ft_selectdata(cfg, dataclean);
    
    % Identify X and B trials in pre that follow A
    pre_trialinfo       = [];
    pre_trialinfo_base  = [];
    previousTrial       = zeros(1,6);
    storyStructure      = unique(data_pre.trialinfo(:,4), 'stable');
    for pretrials = 1:size(storyStructure, 1)
        
        myStory         = storyStructure(pretrials);
        storyIndex      = find(data_pre.trialinfo(:,4) == myStory);
        storyData       = [data_pre.trialinfo(storyIndex,:); previousTrial];
        storyDataPlus   = [previousTrial; data_pre.trialinfo(storyIndex,:)];
        
        findAtrials     = storyDataPlus(:,2) == 1;
        findBtrials     = storyData(:,2) == 2;
        findXtrials     = storyData(:,2) == 3;
        
        validBtrials    = (findAtrials + findBtrials) == 2;
        validXtrials    = (findAtrials + findXtrials) == 2;
        
        % These are all X and B trials that do not follow A
        base_trials     = validBtrials + validXtrials;
        base_pre_trials = base_trials < 1;
        base_pre_trials = storyData(base_pre_trials,:);
        
        pre_trialinfo_base = [pre_trialinfo_base; base_pre_trials];
        
        % Continue with all trials following A
        validBtrials    = storyData(validBtrials,:);
        validXtrials    = storyData(validXtrials,:);
        
        pre_trialinfo   = [pre_trialinfo; validBtrials; validXtrials];
        
        
    end
    
    [~,idx]             = sort(pre_trialinfo(:,1)); 
    pre_trialinfo       = pre_trialinfo(idx,:);
    
    [~,idx]             = sort(pre_trialinfo_base(:,1)); 
    pre_trialinfo_base  = pre_trialinfo_base(idx,:);
    pre_trialinfo_base  = pre_trialinfo_base(13:end,:); %Reject first 12 trials as they are artificial zeros
    
    % Identify X and B trials in post that follow A
    post_trialinfo          = [];
    post_trialinfo_base     = [];
    previousTrial           = zeros(1,6);
    storyStructure          = unique(data_post.trialinfo(:,4), 'stable');
    for posttrials = 1:size(storyStructure, 1)
        
        myStory         = storyStructure(posttrials);
        storyIndex      = find(data_post.trialinfo(:,4) == myStory);
        storyData       = [data_post.trialinfo(storyIndex,:); previousTrial];
        storyDataPlus   = [previousTrial; data_post.trialinfo(storyIndex,:)];
        
        findAtrials     = storyDataPlus(:,2) == 5;
        findBtrials     = storyData(:,2) == 6;
        findXtrials     = storyData(:,2) == 7;
        
        validBtrials    = (findAtrials + findBtrials) == 2;
        validXtrials    = (findAtrials + findXtrials) == 2;
        
        % These are all X and B trials that do not follow A
        base_trials      = validBtrials + validXtrials;
        base_post_trials = base_trials < 1;
        base_post_trials = storyData(base_post_trials,:);
        
        post_trialinfo_base = [post_trialinfo_base; base_post_trials];
        
        % Continue with all trials following A
        validBtrials    = storyData(validBtrials,:);
        validXtrials    = storyData(validXtrials,:);
        
        post_trialinfo  = [post_trialinfo; validBtrials; validXtrials];
        
    end
    
    [~,idx]             = sort(post_trialinfo(:,1)); 
    post_trialinfo      = post_trialinfo(idx,:);
    
    [~,idx]             = sort(post_trialinfo_base(:,1)); 
    post_trialinfo_base = post_trialinfo_base(idx,:);
    post_trialinfo_base = post_trialinfo_base(13:end,:); %Reject first 12 trials as they are artificial zeros
    
    % Identify X trials in pre that are 'surprising'
%     vidcount        = 0;
%     currentTriplet  = [];
%     pre_trialinfo   = [];
%     for pretrials = 1:size(data_pre.trialinfo, 1)
%         
%         currentTriplet  = [currentTriplet; data_pre.trialinfo(pretrials,:), 0];
%         vidcount        = vidcount + 1;
%         
%         if vidcount == 3
%             % If it is a surprisal trial, label it as 1
%             if currentTriplet(2, 3) == 3     
%                 currentTriplet(2,7) = 1;
%                 pre_trialinfo       = [pre_trialinfo; currentTriplet];
%                 vidcount            = 0;
%                 currentTriplet      = [];
%             else
%                 if currentTriplet(3, 3) == 3
%                     currentTriplet(2,7) = 2;
%                 else
%                     currentTriplet(3,7) = 2;
%                 end
%                 pre_trialinfo       = [pre_trialinfo; currentTriplet];
%                 vidcount            = 0;
%                 currentTriplet      = [];
%             end
%         end
%  
%     end
%     
%     % Identify X trials in post that are 'surprising'
%     vidcount         = 0;
%     currentTriplet   = [];
%     post_trialinfo   = [];
%     for posttrials = 1:size(data_post.trialinfo, 1)
%         
%         currentTriplet  = [currentTriplet; data_post.trialinfo(posttrials,:), 0];
%         vidcount        = vidcount + 1;
%         
%         if vidcount == 3
%             % If it is a surprisal trial, label it as 1
%             if currentTriplet(2, 3) == 3     
%                 currentTriplet(2,7) = 1;
%                 post_trialinfo       = [post_trialinfo; currentTriplet];
%                 vidcount            = 0;
%                 currentTriplet      = [];
%             else
%                 if currentTriplet(3, 3) == 3
%                     currentTriplet(2,7) = 2;
%                 else
%                     currentTriplet(3,7) = 2;
%                 end
%                 post_trialinfo       = [post_trialinfo; currentTriplet];
%                 vidcount            = 0;
%                 currentTriplet      = [];
%             end
%         end
%  
%     end
    
    data_pre.trialinfo  = pre_trialinfo;
    data_post.trialinfo = post_trialinfo;
    
    data_pre_base               = data_pre;
    data_post_base              = data_post;
    data_pre_base.trialinfo     = pre_trialinfo_base;
    data_post_base.trialinfo    = post_trialinfo_base;
    
    allTrialsPre    = [allTrialsPre; length(pre_trialinfo)];
    allTrialsPost   = [allTrialsPost; length(post_trialinfo)];
    
    % Select surprisal and expected trials for pre and post that follow A
    cfg             = [];
    %cfg.trials      = find(data_pre.trialinfo(:,7)==2);
    cfg.trials      = find(data_pre.trialinfo(:,2)==2);
    data_pre_exp    = ft_selectdata(cfg, data_pre);
    
    cfg             = [];
    %cfg.trials      = find(data_pre.trialinfo(:,7)==1);
    cfg.trials      = find(data_pre.trialinfo(:,2)==3);
    data_pre_surp   = ft_selectdata(cfg, data_pre);
    
    cfg             = [];
    %cfg.trials      = find(data_post.trialinfo(:,7)==2);
    cfg.trials      = find(data_post.trialinfo(:,2)==6);
    data_post_exp   = ft_selectdata(cfg, data_post);
    
    cfg             = [];
    %cfg.trials      = find(data_post.trialinfo(:,7)==1);
    cfg.trials      = find(data_post.trialinfo(:,2)==7);
    data_post_surp  = ft_selectdata(cfg, data_post);
    
    % Select surprisal and expected trials for pre and post that DO NOT follow A
    cfg                 = [];
    %cfg.trials      = find(data_pre.trialinfo(:,7)==2);
    cfg.trials          = find(data_pre_base.trialinfo(:,2)==2);
    data_pre_exp_base   = ft_selectdata(cfg, data_pre_base);
    
    cfg                 = [];
    %cfg.trials      = find(data_pre.trialinfo(:,7)==1);
    cfg.trials          = find(data_pre_base.trialinfo(:,2)==3);
    data_pre_surp_base  = ft_selectdata(cfg, data_pre_base);
    
    cfg                 = [];
    %cfg.trials      = find(data_post.trialinfo(:,7)==2);
    cfg.trials          = find(data_post_base.trialinfo(:,2)==6);
    data_post_exp_base  = ft_selectdata(cfg, data_post_base);
    
    cfg                 = [];
    %cfg.trials      = find(data_post.trialinfo(:,7)==1);
    cfg.trials          = find(data_post_base.trialinfo(:,2)==7);
    data_post_surp_base = ft_selectdata(cfg, data_post_base);
    
    cfg              = [];
    cfg.feedback     = 'no';
    cfg.method       = 'template';
    cfg.template     = 'CTF275_neighb.mat';
    cfg.planarmethod = 'sincos';
    cfg.channel      = {'MEG'};
    cfg.trials       = 'all';
    cfg.neighbours   = ft_prepare_neighbours(cfg, data_pre_exp);
    
    data_planar_pre_exp     = ft_megplanar(cfg, data_pre_exp);
    data_planar_pre_surp    = ft_megplanar(cfg, data_pre_surp);
    data_planar_post_exp    = ft_megplanar(cfg, data_post_exp);
    data_planar_post_surp   = ft_megplanar(cfg, data_post_surp);
    
    data_planar_pre_exp_base     = ft_megplanar(cfg, data_pre_exp_base);
    data_planar_pre_surp_base    = ft_megplanar(cfg, data_pre_surp_base);
    data_planar_post_exp_base    = ft_megplanar(cfg, data_post_exp_base);
    data_planar_post_surp_base   = ft_megplanar(cfg, data_post_surp_base);
      
    cfg             = [];
    cfg.channel     = 'MEG';
    cfg.trials      = 'all';
    cfg.keeptrials  = 'no';
    erp_pre_exp     = ft_timelockanalysis(cfg, data_planar_pre_exp);
    erp_pre_surp    = ft_timelockanalysis(cfg, data_planar_pre_surp);
    erp_post_exp    = ft_timelockanalysis(cfg, data_planar_post_exp);
    erp_post_surp   = ft_timelockanalysis(cfg, data_planar_post_surp);
    
    erp_pre_exp_base     = ft_timelockanalysis(cfg, data_planar_pre_exp_base);
    erp_pre_surp_base    = ft_timelockanalysis(cfg, data_planar_pre_surp_base);
    erp_post_exp_base    = ft_timelockanalysis(cfg, data_planar_post_exp_base);
    erp_post_surp_base   = ft_timelockanalysis(cfg, data_planar_post_surp_base);
    
    cfg             = [];
    cfg.baseline    = [-0.300 0];
    erp_pre_exp     = ft_timelockbaseline(cfg, erp_pre_exp);
    erp_pre_surp    = ft_timelockbaseline(cfg, erp_pre_surp);
    erp_post_exp    = ft_timelockbaseline(cfg, erp_post_exp);
    erp_post_surp   = ft_timelockbaseline(cfg, erp_post_surp);
    
    erp_pre_exp_base     = ft_timelockbaseline(cfg, erp_pre_exp_base);
    erp_pre_surp_base    = ft_timelockbaseline(cfg, erp_pre_surp_base);
    erp_post_exp_base    = ft_timelockbaseline(cfg, erp_post_exp_base);
    erp_post_surp_base   = ft_timelockbaseline(cfg, erp_post_surp_base);
    
    cfg             = [];
    cfg.method      = 'sum';
    erp_pre_exp     = ft_combineplanar(cfg,erp_pre_exp);
    erp_pre_surp    = ft_combineplanar(cfg,erp_pre_surp);
    erp_post_exp    = ft_combineplanar(cfg,erp_post_exp);
    erp_post_surp   = ft_combineplanar(cfg,erp_post_surp);
    
    erp_pre_exp_base     = ft_combineplanar(cfg,erp_pre_exp_base);
    erp_pre_surp_base    = ft_combineplanar(cfg,erp_pre_surp_base);
    erp_post_exp_base    = ft_combineplanar(cfg,erp_post_exp_base);
    erp_post_surp_base   = ft_combineplanar(cfg,erp_post_surp_base);
    
    all_pre_exp{iSubject}   = erp_pre_exp;
    all_pre_surp{iSubject}  = erp_pre_surp;
    all_post_exp{iSubject}  = erp_post_exp;
    all_post_surp{iSubject} = erp_post_surp;
    
    all_pre_exp_base{iSubject}   = erp_pre_exp_base;
    all_pre_surp_base{iSubject}  = erp_pre_surp_base;
    all_post_exp_base{iSubject}  = erp_post_exp_base;
    all_post_surp_base{iSubject} = erp_post_surp_base;
    
    cfg                     = [];
    cfg.operation           = 'subtract';
    cfg.parameter           = 'avg';
    all_pre_eff{iSubject}   = ft_math(cfg, erp_pre_exp, erp_pre_surp);
    all_post_eff{iSubject}  = ft_math(cfg, erp_post_exp, erp_post_surp);
    all_pre_eff_base{iSubject}   = ft_math(cfg, erp_pre_exp_base, erp_pre_surp_base);
    all_post_eff_base{iSubject}  = ft_math(cfg, erp_post_exp_base, erp_post_surp_base);
    
    all_pre_B_distance{iSubject} = ft_math(cfg, erp_pre_exp_base, erp_pre_exp);
    all_pre_X_distance{iSubject} = ft_math(cfg, erp_pre_surp_base, erp_pre_surp);
    all_post_B_distance{iSubject} = ft_math(cfg, erp_post_exp_base, erp_post_exp);
    all_post_X_distance{iSubject} = ft_math(cfg, erp_post_surp_base, erp_post_surp);
    
    keep all_pre_exp all_pre_surp all_post_exp all_post_surp all_pre_exp_base all_pre_surp_base all_post_exp_base all_post_surp_base subjects root_dir save_dir iSubject mytime all_pre_eff all_post_eff allTrialsPre allTrialsPost all_pre_eff_base all_post_eff_base allTrialsPre_base allTrialsPost_base all_pre_B_distance all_post_B_distance all_pre_X_distance all_post_X_distance
   
end
toc

cfg             = [];
cfg.latency     = [-0.300 2.0];
avg_pre_exp     = ft_timelockgrandaverage(cfg, all_pre_exp{:});
avg_pre_surp    = ft_timelockgrandaverage(cfg, all_pre_surp{:});
avg_post_exp    = ft_timelockgrandaverage(cfg, all_post_exp{:});
avg_post_surp   = ft_timelockgrandaverage(cfg, all_post_surp{:});
avg_pre_eff     = ft_timelockgrandaverage(cfg, all_pre_eff{:});
avg_post_eff    = ft_timelockgrandaverage(cfg, all_post_eff{:});

avg_pre_exp_base     = ft_timelockgrandaverage(cfg, all_pre_exp_base{:});
avg_pre_surp_base    = ft_timelockgrandaverage(cfg, all_pre_surp_base{:});
avg_post_exp_base    = ft_timelockgrandaverage(cfg, all_post_exp_base{:});
avg_post_surp_base   = ft_timelockgrandaverage(cfg, all_post_surp_base{:});
avg_pre_eff_base     = ft_timelockgrandaverage(cfg, all_pre_eff_base{:});
avg_post_eff_base    = ft_timelockgrandaverage(cfg, all_post_eff_base{:});

avg_pre_B_distance  = ft_timelockgrandaverage(cfg, all_pre_B_distance{:});
avg_pre_X_distance  = ft_timelockgrandaverage(cfg, all_pre_X_distance{:});
avg_post_B_distance = ft_timelockgrandaverage(cfg, all_post_B_distance{:});
avg_post_X_distance = ft_timelockgrandaverage(cfg, all_post_X_distance{:});

save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_pre_exp.mat'), 'avg_pre_exp', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_pre_surp.mat'), 'avg_pre_surp', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_post_exp.mat'), 'avg_post_exp', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_post_surp.mat'), 'avg_post_surp', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_pre_eff.mat'), 'avg_pre_eff', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_post_eff.mat'), 'avg_post_eff', '-v7.3')
disp('Done saving 1.')

save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_pre_exp_base.mat'), 'avg_pre_exp_base', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_pre_surp_base.mat'), 'avg_pre_surp_base', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_post_exp_base.mat'), 'avg_post_exp_base', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_post_surp_base.mat'), 'avg_post_surp_base', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_pre_eff_base.mat'), 'avg_pre_eff_base', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_post_eff_base.mat'), 'avg_post_eff_base', '-v7.3')
disp('Done saving 2.')

save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_pre_B_distance.mat'), 'avg_pre_B_distance', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_pre_X_distance.mat'), 'avg_pre_X_distance', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_post_B_distance.mat'), 'avg_post_B_distance', '-v7.3')
save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'avg_post_X_distance.mat'), 'avg_post_X_distance', '-v7.3')
disp('Done saving 3.')

cfg         = [];
cfg.layout  = 'CTF275.lay';
cfg.xlim    = [-0.300 2.0];
figure; ft_multiplotER(cfg, avg_post_surp, avg_post_surp_base);

cfg         = [];
cfg.layout  = 'CTF275.lay';
cfg.xlim    = [0.250 0.600];
figure; ft_topoplotER(cfg, avg_post_eff);


%% ERF Statistics

cfg                         = [];
cfg.method                  = 'template'; 
cfg.template                = 'CTF275_neighb.mat';               
cfg.layout                  = 'CTF275_helmet.mat';                     
cfg.feedback                = 'yes';                            
neighbours                  = ft_prepare_neighbours(cfg, avg_post_exp); 

cfg = [];
cfg.channel                 = {'all'};
cfg.neighbours              = neighbours;
cfg.latency                 = [0 2.000];
cfg.method                  = 'montecarlo';
cfg.statistic               = 'ft_statfun_indepsamplesT';
cfg.correctm                = 'cluster';
cfg.avgovertime             = 'no';
cfg.clusteralpha            = 0.05;
cfg.clusterstatistic        = 'maxsum';
cfg.minnbchan               = 2;
cfg.tail                    = 0;
cfg.clustertail             = 0;
cfg.alpha                   = 0.05;
cfg.numrandomization        = 1000;

subj = length(all_post_exp);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end 
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;


cfg.design = design;
cfg.ivar  = 2;

stat_exp_vs_surp = ft_timelockstatistics(cfg, all_pre_eff{:}, all_post_eff{:});

save(fullfile('/project/3012026.13/processed_RT/surprisal/', 'stat_preeff_vs_posteff.mat'), 'stat_exp_vs_surp', '-v7.3')
disp('Done saving.')    


%% Plot cluster

cfg = [];
cfg.alpha = 0.05;
ft_clusterplot(cfg, stat_exp_vs_surp)
ft_hastoolbox('brewermap', 1);
set(gcf,'color','w')
colormap(flipud(brewermap(64,'RdBu'))) 



indexProb = stat_exp_vs_surp.prob <= 0.05;
indexProb = sum(indexProb,2);
indexProb = indexProb > 0;
probLabels = stat_exp_vs_surp.label(indexProb);


% Plot channels under specifified p-value
cfg         = [];
cfg.layout  = 'CTF275.lay';
cfg.xlim    = [1.3613 1.4046];
cfg.highlight = 'on';
cfg.highlightchannel = probLabels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 30;
cfg.comment = 'no';
cfg.style = 'blank';
figure; ft_topoplotER(cfg, avg_post_eff);
ft_hastoolbox('brewermap', 1);
set(gcf,'color','w')
colormap(flipud(brewermap(64,'RdBu'))) 

%% Barplot


index_time1 = 1.3613;
[val,idx]=min(abs(avg_pre_exp.time-index_time1));
index_time1= idx;

index_time2 = 1.4046;
[val,idx]=min(abs(avg_pre_exp.time-index_time2));
index_time2= idx;

mean_pre_exp    = mean(mean(avg_pre_exp.avg(indexProb,index_time1:index_time2),1));
mean_pre_surp   = mean(mean(avg_pre_surp.avg(indexProb,index_time1:index_time2),1));
mean_post_exp   = mean(mean(avg_post_exp.avg(indexProb,index_time1:index_time2),1));
mean_post_surp  = mean(mean(avg_post_surp.avg(indexProb,index_time1:index_time2),1));

std_pre_exp    = mean(mean(sqrt(avg_pre_exp.var(indexProb,index_time1:index_time2)),1));
std_pre_surp   = mean(mean(sqrt(avg_pre_surp.var(indexProb,index_time1:index_time2)),1));
std_post_exp   = mean(mean(sqrt(avg_post_exp.var(indexProb,index_time1:index_time2)),1));
std_post_surp  = mean(mean(sqrt(avg_post_surp.var(indexProb,index_time1:index_time2)),1));

mean_pre_B_distance     = mean(mean(avg_pre_B_distance.avg(indexProb,index_time1:index_time2),1));
mean_pre_X_distance     = mean(mean(avg_pre_X_distance.avg(indexProb,index_time1:index_time2),1));
mean_post_B_distance    = mean(mean(avg_post_B_distance.avg(indexProb,index_time1:index_time2),1));
mean_post_X_distance    = mean(mean(avg_post_X_distance.avg(indexProb,index_time1:index_time2),1));

std_pre_B_distance      = mean(mean(sqrt(avg_pre_B_distance.var(indexProb,index_time1:index_time2)),1));
std_pre_X_distance      = mean(mean(sqrt(avg_pre_X_distance.var(indexProb,index_time1:index_time2)),1));
std_post_B_distance     = mean(mean(sqrt(avg_post_B_distance.var(indexProb,index_time1:index_time2)),1));
std_post_X_distance     = mean(mean(sqrt(avg_post_X_distance.var(indexProb,index_time1:index_time2)),1));

mean_pre_exp_base    = mean(mean(avg_pre_exp_base.avg(indexProb,index_time1:index_time2),1));
mean_pre_surp_base   = mean(mean(avg_pre_surp_base.avg(indexProb,index_time1:index_time2),1));
mean_post_exp_base   = mean(mean(avg_post_exp_base.avg(indexProb,index_time1:index_time2),1));
mean_post_surp_base  = mean(mean(avg_post_surp_base.avg(indexProb,index_time1:index_time2),1));

std_pre_exp_base    = mean(mean(sqrt(avg_pre_exp_base.var(indexProb,index_time1:index_time2)),1));
std_pre_surp_base   = mean(mean(sqrt(avg_pre_surp_base.var(indexProb,index_time1:index_time2)),1));
std_post_exp_base   = mean(mean(sqrt(avg_post_exp_base.var(indexProb,index_time1:index_time2)),1));
std_post_surp_base  = mean(mean(sqrt(avg_post_surp_base.var(indexProb,index_time1:index_time2)),1));


%% Bar plot 1

X = categorical({'Pre B_A','Pre X_A','Post B_A','Post X_A','Pre B_X','Pre X_B','Post B_X','Post X_B'});
X = reordercats(X,{'Pre B_A','Pre X_A','Post B_A','Post X_A','Pre B_X','Pre X_B','Post B_X','Post X_B'});
mean_amplitude = [mean_pre_exp, mean_pre_surp, mean_post_exp, mean_post_surp, mean_pre_exp_base, mean_pre_surp_base, mean_post_exp_base, mean_post_surp_base];
mean_std = [std_pre_exp, std_pre_surp, std_post_exp, std_post_surp, std_pre_exp_base, std_pre_surp_base, std_post_exp_base, std_post_surp_base];

hold on
hb = bar(X, mean_amplitude, 'FaceColor','flat')
errorbar(X,mean_amplitude,mean_std,'.')
hb.CData(5,:) = [0 0.8 0.8];
hb.CData(6,:) = [0 0.8 0.8];
hb.CData(7,:) = [0 0.8 0.8];
hb.CData(8,:) = [0 0.8 0.8];
title('Mean Amplitude for 1.3613s - 1.4046s')
ylabel('Planar gradient (fT/cm)')
set(gcf,'color','w')


%% Bar plot 2

X = categorical({'Pre (BX-BA)','Post (BX-BA)','Pre (XB-XA)','Post (XB-XA)'});
X = reordercats(X,{'Pre (BX-BA)','Post (BX-BA)','Pre (XB-XA)','Post (XB-XA)'});
mean_amplitude = [mean_pre_B_distance, mean_post_B_distance, mean_pre_X_distance, mean_post_X_distance];
mean_std = [std_pre_B_distance, std_post_B_distance, std_pre_X_distance, std_post_X_distance];

hold on
hb = bar(X, mean_amplitude, 'FaceColor','flat')
errorbar(X,mean_amplitude,mean_std,'.')
hb.CData(2,:) = [0 0.8 0.8];
hb.CData(4,:) = [0 0.8 0.8];
title('Mean Amplitude for 1.3613s - 1.4046s')
ylabel('Planar gradient (fT/cm)')
set(gcf,'color','w')

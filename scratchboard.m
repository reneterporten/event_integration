% Average over subjects
cfg1 = [];
cfg1.appenddim = 'rpt';

cfg2 = [];
cfg2.avgoverrpt = 'yes';

for m = 1:size(T,1)
    %Favg{m,1} = ft_selectdata(cfg2, ft_appendtimelock(cfg1, F{m,:}));
    %Favg{m,1} = ft_selectdata(cfg2, ft_appendfreq(cfg1, F{m,:}));
    Tavg{m,1} = ft_selectdata(cfg2, ft_appendtimelock(cfg1, T{m,:}));
end

% ensure same time axis (subsampling misalignment may occur due to
% resampling step
for k = 1:numel(Tavg)
  %Favg{k,1}.time = Favg{1}.time;
  Tavg{k,1}.time = Tavg{1}.time;
end

cfg1 = [];
cfg1.appenddim = 'rpt';

cfg2 = [];
cfg2.avgoverrpt = 'yes';

% Average over AB pre
FABpre = ft_selectdata(cfg2, ft_appendtimelock(cfg1, Favg{1}, Favg{2}));
%FABpre = ft_selectdata(cfg2, ft_appendfreq(cfg1, Favg{1}, Favg{2}));

% Average over AB post
FABpost = ft_selectdata(cfg2, ft_appendtimelock(cfg1, Favg{4}, Favg{5}));
%FABpost = ft_selectdata(cfg2, ft_appendfreq(cfg1, Favg{4}, Favg{5}));

Fdiff = FABpost;
%Fdiff.powspctrm = Favg{6}.powspctrm - FABpost.powspctrm;
%Fdiff.trial = Favg{6}.trial - FABpost.trial;
Fdiff.trial = Favg{6}.trial - Favg{5}.trial;

cfg = [];
cfg.viewmode = 'butterfly';
cfg.allowoverlap = 'yes';
ft_databrowser(cfg, Fdiff)
set(gcf,'color','w');

cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.xlim = [0 0.5];
cfg.ylim = [-3.6e-15 2.7e-15];
cfg.comment = 'no';
ft_topoplotER(cfg, Fdiff);
set(gcf,'color','w');
cfg.xlim = [0.5 1.0];
cfg.ylim = [-3.6e-15 2.7e-15];
cfg.comment = 'no';
ft_topoplotER(cfg, Fdiff);
set(gcf,'color','w');
cfg.xlim = [1.0 1.5];
cfg.ylim = [-3.6e-15 2.7e-15];
cfg.comment = 'no';
ft_topoplotER(cfg, Fdiff);
set(gcf,'color','w');
cfg.xlim = [1.5 2.0];
cfg.ylim = [-3.6e-15 2.7e-15];
cfg.comment = 'no';
ft_topoplotER(cfg, Fdiff);
set(gcf,'color','w');

cfg = [];
cfg.layout = 'CTF275_helmet.mat';
%ft_multiplotER(cfg, Tavg{20}, Tavg{8})
cfg.linewidth = 1.0;
%ft_multiplotER(cfg, Tavg{1}, Tavg{7}, Tavg{12},Tavg{16}, Tavg{19}, Tavg{21})
ft_multiplotER(cfg, Tavg{7}, Tavg{19})
set(gcf,'color','w');

% find sign channels
signchan= stat_BXpost_BXpre.negclusterslabelmat(:,73:139); %sign timewindow
signchan = sum(signchan, 2);
signchan = signchan(:) > 10;
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
%cfg.channel = Tavg{1,1}.label(signchan);
fg.channel = {'MLO13', 'MLO14', 'MLO24', 'MLP23', 'MLP33', 'MLP34', 'MLP35', 'MLP42', 'MLP43', 'MLP44', 'MLP45', 'MLP53', 'MLP54', 'MLP55', 'MLP56', 'MLP57', 'MLT14', 'MLT15', 'MLT16', 'MLT25', 'MLT26', 'MLT27', 'MLT36', 'MLT37', 'MLT46'};
ft_singleplotER(cfg, Tavg{20}, Tavg{8})
set(gcf,'color','w');

cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.xlim = [0 0.5];
cfg.comment = 'no';
cfg.highlight = 'on';
cfg.highlightsymbol = '.';
cfg.highlightsize = 24;
cfg.highlightchannel = Tavg{1,1}.label(signchan);
ft_topoplotER(cfg, Tavg{20});



cfg = [];
cfg.layout = 'CTF275_helmet.mat';
%cfg.xlim = [0.45 1.1];
figure
ft_singleplotTFR(cfg, Fdiff);


%% Check trial order across subjects


%ignSub = [4, 10, 12, 20, 30, 32];
ignSub = [999];
countSubs = 1;
allpos = [];
for iSubject = 1:length(subjects)
    
    if ~ismember(iSubject, ignSub)
        
        disp(subjects{iSubject})
        load(fullfile(root_dir, [subjects{iSubject},'_dataclean.mat']), 'dataclean');

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
        
        datatrials = dataclean.trialinfo(dataclean.trialinfo(:,2) ~= 4,:);
        
        % Check overall position
        % First position
        first_vec = datatrials(1:3:end,2);
        % Second position
        second_vec = datatrials(2:3:end,2);
        % Third position
        third_vec = datatrials(3:3:end,2);
        
        Apos = [sum(sum(first_vec == [1, 5],2)), sum(sum(second_vec == [1, 5],2)), sum(sum(third_vec == [1, 5],2))];
        Bpos = [sum(sum(first_vec == [2, 6],2)), sum(sum(second_vec == [2, 6],2)), sum(sum(third_vec == [2, 6],2))];
        Xpos = [sum(sum(first_vec == [3, 7],2)), sum(sum(second_vec == [3, 7],2)), sum(sum(third_vec == [3, 7],2))];
        
        Apospre = [sum(sum(first_vec == [1],2)), sum(sum(second_vec == [1],2)), sum(sum(third_vec == [1],2))];
        Bpospre = [sum(sum(first_vec == [2],2)), sum(sum(second_vec == [2],2)), sum(sum(third_vec == [2],2))];
        Xpospre = [sum(sum(first_vec == [3],2)), sum(sum(second_vec == [3],2)), sum(sum(third_vec == [3],2))];
        
        Apospost = [sum(sum(first_vec == [5],2)), sum(sum(second_vec == [5],2)), sum(sum(third_vec == [5],2))];
        Bpospost = [sum(sum(first_vec == [6],2)), sum(sum(second_vec == [6],2)), sum(sum(third_vec == [6],2))];
        Xpospost = [sum(sum(first_vec == [7],2)), sum(sum(second_vec == [7],2)), sum(sum(third_vec == [7],2))];
        
        allpos.pos(countSubs,:) = [Apos, Bpos, Xpos];
        allpos.pospre(countSubs,:) = [Apospre, Bpospre, Xpospre];
        allpos.pospost(countSubs,:) = [Apospost, Bpospost, Xpospost];
        allpos.subs{countSubs} = subjects{iSubject};
        allpos.posord = {'A1', 'A2', 'A3', 'B1', 'B2', 'B3','X1', 'X2', 'X3'};
        
        countSubs = countSubs +1;
    end

end

% Line plot
figure;
positions = [1, 4, 7];
for i = 1:size(allpos.pos,1)
    for x = 1:3
        subplot(3,1,x)
        hold on
        plot(allpos.pos(i,positions(x):positions(x)+2))
        xticks(1:3)
        xticklabels(allpos.posord(positions(x):positions(x)+2))
    end
end

% Image plot
figure; hAxes = gca;
clims = [28 68];
imagesc(hAxes, allpos.pos, clims)
xticks(1:9)
xticklabels(allpos.posord)
yticks(1:length(allpos.pos(:,1)))
yticklabels(allpos.subs)
ft_hastoolbox('brewermap', 1);         
colormap(hAxes, flipud(brewermap(64,'RdBu')))
set(gcf,'color','w');
title('Position Count for all Trials')
colorbar

% Pre
figure; hAxes = gca;
clims = [14 34];
imagesc(hAxes, allpos.pospre, clims)
xticks(1:9)
xticklabels(allpos.posord)
yticks(1:length(allpos.pospre(:,1)))
yticklabels(allpos.subs)
ft_hastoolbox('brewermap', 1);         
colormap(hAxes, flipud(brewermap(64,'RdBu')))
set(gcf,'color','w');
title('Position Count for Pre Trials')
colorbar

% Post
figure; hAxes = gca;
clims = [14 34];
imagesc(hAxes, allpos.pospost, clims)
xticks(1:9)
xticklabels(allpos.posord)
yticks(1:length(allpos.pospost(:,1)))
yticklabels(allpos.subs)
ft_hastoolbox('brewermap', 1);         
colormap(hAxes, flipud(brewermap(64,'RdBu')))
set(gcf,'color','w');
title('Position Count for Post Trials')
colorbar

% Subject average
% All Trials
avgpos = mean(allpos.pos,1);
figure;
bar(avgpos)
xticks(1:9)
xticklabels(allpos.posord)
ylim([0 60])
title('Position Count-Average for all Trials')

% Pre
avgpospre = mean(allpos.pospre,1);
figure;
bar(avgpospre)
xticks(1:9)
xticklabels(allpos.posord)
ylim([0 30])
title('Position Count-Average for Pre Trials')

% Post
avgpospost = mean(allpos.pospost,1);
figure;
bar(avgpospost)
xticks(1:9)
xticklabels(allpos.posord)
ylim([0 30])
title('Position Count-Average for Post Trials')


%%
ignSub = [4, 10, 12, 20, 30, 32];
for m = 1:numel(subjects)
    subj = subjects{m}; 
    if ~ismember(m, ignSub)
        m
        qsubfeval('rt_sensorlevelanalysis', subj, 'saveflag', '/project/3012026.13/jansch/', 'type', 'timelock123', 'memreq', (1024^3)*12, 'timreq', 60*60)
    end
end


%%
suff = '123timelock123.mat';
datadir = '/project/3012026.13/jansch';
cd(datadir);

cfg1 = [];
cfg1.appenddim = 'rpt';

cfg2 = [];
cfg2.avgoverrpt = 'yes';
d = dir(sprintf('sub*%s*',suff));

for k = 1:numel(d)
      
  disp(d(k).name)
  tmp = load(d(k).name);
  varname = fieldnames(tmp);
  data = tmp.(varname{1});
  
  istimelock = ft_datatype(data(1), 'timelock');
  
  for m = 1:numel(data)
    
    % these are the indivdual pairwise correlations
    F{m,k} = data(m);
    F{m,k}.time = F{1,1}.time;
  end

end


%%

cfg                 = [];
cfg.xlim            = [0.6 0.8];
%cfg.zlim            = [-0.005 0.005];
cfg.layout          = 'CTF275_helmet.mat';
cfg.parameter       = 'stat';
cfg.comment         = 'xlim';
cfg.commentpos      = 'title';
cfg.marker          = 'off';
cfg.colorbar        = 'no';

figure; ft_topoplotER(cfg, stat_23); colorbar
set(gcf,'color','w')
ft_hastoolbox('brewermap', 1);         
colormap(flipud(brewermap(64,'RdBu')))


%%

ignSub = [4, 10, 12, 20, 30, 32];

load(fullfile(root_dir, [subjects{10},'_dataclean.mat']), 'dataclean');

trialidx = ones(length(dataclean.trialinfo),1);
xend = 42;
xstart = 1;

storvec = ones(12,2);
for stor = 1:12
    %storysel = xstart:(stor*xend);
    %xstart = xstart + xend;
    storvec(stor,:) = [stor, length(find(dataclean.trialinfo(:,4) == stor))];
end

% Subjects easily fixable:
% 4, trl 388
% 10, trl 107
% 12, trl 296
% 20, trl ??
% 30, trl 296
% 32, trl 103


%% Some subjects have an additional anomalous trial that needs to be rejected

probSubs = {'sub-007'; 'sub-013'; 'sub-016'; 'sub-036'; 'sub-038'};
probTrials = [388, 107, 296, 296, 103];

if any(strcmp(probSubs, subj))
    subIdx = find(strcmp(probSubs, subj));
    cfg = [];
    cfg.trials = find(dataclean.trialinfo(:,1) ~= probTrials(subIdx));
    dataclean = ft_selectdata(cfg, dataclean);
end






















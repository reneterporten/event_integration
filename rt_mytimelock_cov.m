function [tlck] = rt_mytimelock_cov(root_dir, subj)

load(fullfile(root_dir, [subj,'_dataclean.mat']));

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

cfg             = [];
cfg.demean      = 'yes';
dataclean       = ft_preprocessing(cfg, dataclean);

cfg             = [];
cfg.covariance  = 'yes';
tlck            = ft_timelockanalysis(cfg, dataclean);

cfg             = [];
cfg.latency     = [0 2.0];
tlck            = ft_selectdata(cfg,tlck);

save(fullfile('/project/3012026.13/processed_RT/source reconstruction/', subj, 'tlck_filter.mat'), 'tlck')

function TF_basecont_qsubRT(save_dir, subj, data_A_pre, data_A_post, data_B_pre, data_B_post, data_X_pre, data_X_post)

cfg             = [];
cfg.output      = 'pow';
cfg.channel     = 'MEG';
cfg.method      = 'mtmconvol';
cfg.taper       = 'hanning';
cfg.foi         = 2:1:40;
cfg.trials      = 'all';
cfg.t_ftimwin   = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec % 1.0 sec
cfg.toi         = [-1:0.05:2];
cfg.pad         ='nextpow2';
cfg.keeptrials  = 'yes'; %% we keep the single trials
        
freq_A_pre      = ft_freqanalysis(cfg, data_A_pre);
freq_A_post     = ft_freqanalysis(cfg, data_A_post);
freq_B_pre      = ft_freqanalysis(cfg, data_B_pre);
freq_B_post     = ft_freqanalysis(cfg, data_B_post);
freq_X_pre      = ft_freqanalysis(cfg, data_X_pre);
freq_X_post     = ft_freqanalysis(cfg, data_X_post);

cfg             = [];
cfg.method      = 'sum';
freq_A_pre      = ft_combineplanar(cfg,freq_A_pre);
freq_A_post     = ft_combineplanar(cfg,freq_A_post);
freq_B_pre      = ft_combineplanar(cfg,freq_B_pre);
freq_B_post     = ft_combineplanar(cfg,freq_B_post);
freq_X_pre      = ft_combineplanar(cfg,freq_X_pre);
freq_X_post     = ft_combineplanar(cfg,freq_X_post);

save(fullfile(save_dir, [subj,'_TF_basecon']),'freq_A_pre','freq_A_post','freq_B_pre','freq_B_post', 'freq_X_pre','freq_X_post')

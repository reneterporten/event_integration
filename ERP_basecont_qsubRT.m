function ERP_basecont_qsubRT(save_dir, subj, data_A_pre, data_A_post, data_B_pre, data_B_post, data_X_pre, data_X_post)

cfg             = [];
cfg.channel     = 'MEG';
cfg.trials      = 'all';
cfg.keeptrials  = 'yes';
erp_A_pre      = ft_timelockanalysis(cfg, data_A_pre);
erp_A_post     = ft_timelockanalysis(cfg, data_A_post);
erp_B_pre      = ft_timelockanalysis(cfg, data_B_pre);
erp_B_post     = ft_timelockanalysis(cfg, data_B_post);
erp_X_pre      = ft_timelockanalysis(cfg, data_X_pre);
erp_X_post     = ft_timelockanalysis(cfg, data_X_post);

cfg             = [];
cfg.method      = 'sum';
erp_A_pre      = ft_combineplanar(cfg,erp_A_pre);
erp_A_post     = ft_combineplanar(cfg,erp_A_post);
erp_B_pre      = ft_combineplanar(cfg,erp_B_pre);
erp_B_post     = ft_combineplanar(cfg,erp_B_post);
erp_X_pre      = ft_combineplanar(cfg,erp_X_pre);
erp_X_post     = ft_combineplanar(cfg,erp_X_post);

cfg             = [];
cfg.latency     = [0 2.0];
erp_A_pre      = ft_selectdata(cfg,erp_A_pre);
erp_A_post     = ft_selectdata(cfg,erp_A_post);
erp_B_pre      = ft_selectdata(cfg,erp_B_pre);
erp_B_post     = ft_selectdata(cfg,erp_B_post);
erp_X_pre      = ft_selectdata(cfg,erp_X_pre);
erp_X_post     = ft_selectdata(cfg,erp_X_post);

save(fullfile(save_dir, [subj,'_ERP_basecon']),'erp_A_pre','erp_A_post','erp_B_pre','erp_B_post', 'erp_X_pre','erp_X_post')
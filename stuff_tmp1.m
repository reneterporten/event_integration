%% 
% get the data
%cfg      = ft_getopt(varargin, 'cfg_preproc', []);

subj = 'sub-004';
cfg = [];
root_dir = '/project/3012026.13/processed';
data     = rt_mytimelockv3(root_dir, subj, 'cfg_preproc', cfg);


data = ft_appenddata([], data{:});

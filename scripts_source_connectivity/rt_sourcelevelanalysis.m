function rt_sourcelevelanalysis(subj, varargin)

% Function to define leadfield and estimate activity at source level
% Returns connectivity measures between the hippocampus and the mPFC

if nargin<1 || isempty(subj)
  subj = 'sub-004';
end

saveflag    = ft_getopt(varargin, 'saveflag', false);
savepath    = ft_getopt(varargin, 'savepath', '/project/3012026.13/jansch/');
savename    = ft_getopt(varargin, 'savename', 'coh'); 
headdata    = ft_getopt(varargin, 'headdata', fullfile('/project/3012026.13/jansch/', strcat(subj, '_headsource.mat')));
cfgfreq     = ft_getopt(varargin, 'cfgfreq', []);
comp        = ft_getopt(varargin, 'comp', 'abxprepost');
load(headdata) % Loads headmodel and segmented mri


%% Prepare the sensor level data for source analysis

root_dir = '/project/3012026.13/processed';
cfgpreproc.lpfilter = 'no';
data     = rt_mytimelockv3(root_dir, subj, 'cfg_preproc', cfgpreproc, 'latency', [0.3 2], 'doplanar', 0);
data     = ft_appenddata([], data{:});
data     = removefields(data, {'elec' 'cfg'});


%% Do spectral transformation for the theta band

cfg             = [];
cfg.method      = 'mtmfft';
cfg.output      = 'fourier';
cfg.foilim      = ft_getopt(cfgfreq, 'foilim', [5.5 5.5]);
cfg.tapsmofrq   = ft_getopt(cfgfreq, 'tapsmofrq', 1.5);
cfg.pad         = 4;
cfg.trials      = 'all';
freq            = ft_freqanalysis(cfg, data);


%% Prepare leadfield
cfg                     = [];
cfg.grad                = freq.grad;
cfg.headmodel           = headmodel;
cfg.sourcemodel         = sourcemodel;
cfg.channel             = {'MEG'};
cfg.singleshell.batchsize = 2000;
sourcemodel             = ft_prepare_leadfield(cfg);


%% Source analysis

freqcsd = ft_checkdata(freq, 'cmbstyle', 'fullfast');

cfg                     = [];
cfg.method              = 'dics';
cfg.frequency           = freqcsd.freq;
cfg.sourcemodel         = sourcemodel;
cfg.headmodel           = headmodel;
cfg.dics.projectnoise   = 'yes';
cfg.dics.lambda         = '20%';
cfg.dics.keepfilter     = 'yes';
cfg.dics.fixedori       = 'yes';
cfg.dics.weightnorm     = 'unitnoisegain';

source = ft_sourceanalysis(cfg, freqcsd);


%% Apply atlas information to source level data

atlasfile = fullfile('/project/3012026.13/jansch/brainnetome_atlas_grid');
load(atlasfile);

atlas_grid.pos = source.pos;

cfg                 = [];
cfg.parcellation    = 'tissue';
cfg.method          = 'svd';
cfg.numcomponent    = 5;
% Create virtual channels from conditions
switch comp
    case 'prepost'  
        conditions  = [1 2]; %pre post
        conlabel    = {'pre', 'post'};
        for c = 1:numel(conditions)
            % Select data from freq and create virtual channels
            cfgsel          = [];
            cfgsel.trials   = find(freq.trialinfo(:,5)==conditions(c));
            freqsel         = ft_selectdata(cfgsel, freq);
            vc{c}           = ft_virtualchannel(cfg, freqsel, source, atlas_grid);
        end   
    case 'abxprepost'
        conditions  = [1 2 3 5 6 7]; % a b x (pre), a b x (post)
        conlabel    = {'a_pre', 'b_pre', 'x_pre', 'a_post', 'b_post', 'x_post'};
        for c = 1:numel(conditions)
            % Select data from freq and create virtual channels
            cfgsel          = [];
            cfgsel.trials   = find(freq.trialinfo(:,2)==conditions(c));
            freqsel         = ft_selectdata(cfgsel, freq);
            vc{c}           = ft_virtualchannel(cfg, freqsel, source, atlas_grid);
        end 
end


%% Connectivity analysis (coherence)

for k = 1:numel(vc)
    vccsd       = ft_checkdata(vc{k}, 'cmbstyle', 'fullfast');

    cfg         = [];
    cfg.method  = 'coh'; % for instance
    cfg.complex = 'abs';
    coh{k}         = ft_connectivityanalysis(cfg, vccsd);

    % get the index grouping of the ROI's components
    label = coh{k}.label;
    for j = 1:numel(label)
      label{j} = label{j}(1:end-4);
    end
    [ulabel, i1, i2] = unique(label, 'stable');

    % create a projection matrix for fast averaging
    P = sparse(i2, (1:numel(label))', ones(numel(label),1));
    P = P./sum(P,2);

    cfg             = [];
    cfg.method      = 'mim';
    cfg.indices     = i2(:);
    mim{k}          = ft_connectivityanalysis(cfg, vccsd);

    cfg             = [];
    cfg.method      = 'coh';
    cfg.complex     = 'absimag';
    imcoh{k}        = ft_connectivityanalysis(cfg, vccsd);

    coh{k}.cohspctrm   = P*coh{k}.cohspctrm*P';
    imcoh{k}.cohspctrm = P*imcoh{k}.cohspctrm*P';
    coh{k}.label       = ulabel;
    imcoh{k}.label     = ulabel;
    mim{k}.label       = ulabel; % for readability: check whether this is correct, i.e. that it doesn't mix up the labels
end


%% Save variables

if saveflag
    fname = fullfile(savepath, sprintf('%s_%s', subj, savename));
    save(fname, 'coh', 'imcoh', 'mim', 'conlabel');
end



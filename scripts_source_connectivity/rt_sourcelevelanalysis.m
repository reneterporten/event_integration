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
% Per condition
cfg.trials      = find(data.trialinfo(:,5) == 1); %pre
freqpre         = ft_freqanalysis(cfg, data);
cfg.trials      = find(data.trialinfo(:,5) == 2); %post
freqpost        = ft_freqanalysis(cfg, data);


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
vc                  = ft_virtualchannel(cfg, freq, source, atlas_grid);
vcpre               = ft_virtualchannel(cfg, freqpre, source, atlas_grid);
vcpost              = ft_virtualchannel(cfg, freqpost, source, atlas_grid);


%% Connectivity analysis (coherence)

vccsd       = ft_checkdata(vc, 'cmbstyle', 'fullfast');
vcprecsd    = ft_checkdata(vcpre, 'cmbstyle', 'fullfast');
vcpostcsd   = ft_checkdata(vcpost, 'cmbstyle', 'fullfast');

cfg         = [];
cfg.method  = 'coh'; % for instance
cfg.complex = 'abs';
coh         = ft_connectivityanalysis(cfg, vccsd); % pre
cohpre      = ft_connectivityanalysis(cfg, vcprecsd); % pre
cohpost     = ft_connectivityanalysis(cfg, vcpostcsd); % post

% get the index grouping of the ROI's components
label = coh.label;
for k = 1:numel(label)
  label{k} = label{k}(1:end-4);
end
[ulabel, i1, i2] = unique(label, 'stable');

% create a projection matrix for fast averaging
P = sparse(i2, (1:numel(label))', ones(numel(label),1));
P = P./sum(P,2);

cfg         = [];
cfg.method  = 'mim';
cfg.indices = i2(:);
mim         = ft_connectivityanalysis(cfg, vccsd);

cfg         = [];
cfg.method  = 'coh';
cfg.complex = 'absimag';
imcoh       = ft_connectivityanalysis(cfg, vccsd);

coh.cohspctrm   = P*coh.cohspctrm*P';
imcoh.cohspctrm = P*imcoh.cohspctrm*P';
coh.label       = ulabel;
imcoh.label     = ulabel;
mim.label       = ulabel; % for readability: check whether this is correct, i.e. that it doesn't mix up the labels


%% Save variables

if saveflag
    fname = fullfile(savepath, sprintf('%s_%s', subj, savename));
    save(fname, 'cohpre', 'cohpost');
end



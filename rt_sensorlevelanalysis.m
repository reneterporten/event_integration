function varargout = rt_sensorlevelanalysis(subj, cfg_preproc, type, mccaflag, saveflag)

if nargin<1 || isempty(subj)
  subj = 'sub-004';
end

if nargin<2 || isempty(cfg_preproc)
  cfg = [];
else
  cfg = cfg_preproc;
end

if nargin<3 || isempty(type)
  type = 'timelock';
end

if nargin<4 || isempty(mccaflag)
  mccaflag = false;
end

if nargin<5 || isempty(saveflag)
  saveflag = false;
else
  if ischar(saveflag)
    savepath = saveflag;
    saveflag = true;
  end
end

if saveflag
  if ~exist('savepath', 'var')
    savepath = '';
  end
end

if isempty(cfg) && ~strcmp(type, 'tfr')
  cfg.lpfilter   = 'yes';
  cfg.lpfreq     = 35;
  cfg.lpfilttype = 'firws';
elseif isempty(cfg) && strcmp(type, 'tfr')
  cfg.lpfilter = 'no';
  latency      = [-1 2-1./300];
end

if ~exist('latency', 'var')
  latency = [];
end

root_dir = '/project/3012026.13/processed';
data     = rt_mytimelockv3(root_dir, subj, 'cfg_preproc', cfg, 'latency', latency);
data     = ft_appenddata([], data{:});
data     = removefields(data, {'elec' 'cfg'});

if mccaflag
  data = rt_mcca(data);
end

switch type
  case 'timelock'
    conds = [1 2 3 5 6 7];
    for k = 1:numel(conds)
      tmpcfg        = [];
      tmpcfg.trials = find(data.trialinfo(:,2)==conds(k));
      tmpcfg.preproc.baselinewindow = [-0.1 0];
      tmpcfg.preproc.demean         = 'yes';
      tlck(k) = ft_timelockanalysis(tmpcfg, data);
    end
    
    for k = 1:6
      tlck(k).time = tlck(1).time;
    end
    
    if ~mccaflag
      % combine planar gradients
      for k = 1:numel(tlck)
        tlck_out(k) = ft_combineplanar([], tlck(k));
      end
      tlck = tlck_out;
      clear tlck_out;
    end
    
    varargout{1} = tlck;
    if saveflag
      if mccaflag
        suff = '_mcca';
      else
        suff = '';
      end
      tlck = rmfield(tlck, 'cfg'); % this one is really big due to rt_mytimelockv3
      
      fname = fullfile(savepath, sprintf('%s_tlck%s', subj, suff));
      save(fname, 'tlck');
      clear data
    end
  case 'tfr'
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.foi    = 2:2:30;
    cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.5;
    cfg.taper  = 'hanning';
    cfg.toi    = (-150:15:450)./300;
    cfg.pad    = 4;
    
    conds = [1 2 3 5 6 7];
    for k = 1:numel(conds)
      cfg.trials = find(data.trialinfo(:,2)==conds(k));
      if ~mccaflag
        freq(k)    = ft_combineplanar([], ft_freqanalysis(cfg, data));
      else
        freq(k)    = ft_freqanalysis(cfg, data);
      end
    end
    
    varargout{1} = freq;
    if saveflag
      if mccaflag
        suff = '_mcca';
      else
        suff = '';
      end
      %freq = rmfield(freq, 'cfg'); % this one is really big due to rt_mytimelockv3
      
      fname = fullfile(savepath, sprintf('%s_freq%s', subj, suff));
      save(fname, 'freq');
      clear data
    end
    
    
  case 'rsa'
    % to-be-implemented
end

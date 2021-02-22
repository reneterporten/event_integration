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

if isempty(cfg)
  cfg.lpfilter   = 'yes';
  cfg.lpfreq     = 35;
  cfg.lpfilttype = 'firws';
end

root_dir = '/project/3012026.13/processed';
data     = rt_mytimelockv3(root_dir, subj, 'cfg_preproc', cfg);
data     = ft_appenddata([], data{:});

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
    
    varargout{1} = tlck;
    if saveflag
      if mccaflag
        suff = '_mcca';
      else
        suff = '';
      end
      fname = fullfile(savepath, sprintf('%s_tlck%s', subj, suff));
      save(fname, 'tlck');
    end
  case 'tfr'
    % to-be-implemented
  case 'rsa'
    % to-be-implemented
end

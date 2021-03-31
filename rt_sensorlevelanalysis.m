function varargout = rt_sensorlevelanalysis(subj, varargin)

if nargin<1 || isempty(subj)
  subj = 'sub-004';
end

cfgpreproc = ft_getopt(varargin, 'cfgpreproc', []);
cfgfreq    = ft_getopt(varargin, 'cfgfreq', []);
type       = ft_getopt(varargin, 'type', 'timelock');
mccaflag   = ft_getopt(varargin, 'mccaflag', false);
saveflag   = ft_getopt(varargin, 'saveflag', false); % can be string to path, or boolean
savename   = ft_getopt(varargin, 'savename', '');
poscorrect = ft_getopt(varargin, 'poscorrect', false);

if isempty(saveflag)
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
  if isempty(savename)
    savename = type;
  end
end

if isempty(cfgpreproc) && ~startsWith(type, 'tfr')
  cfgpreproc.lpfilter   = 'yes';
  cfgpreproc.lpfreq     = 35;
  cfgpreproc.lpfilttype = 'firws';
elseif isempty(cfgpreproc) && startsWith(type, 'tfr')
  cfgpreproc.lpfilter = 'no';
  latency      = [-1 2-1./300];
end

if ~exist('latency', 'var')
  latency = [];
end

root_dir = '/project/3012026.13/processed';
data     = rt_mytimelockv3(root_dir, subj, 'cfg_preproc', cfgpreproc, 'latency', latency);
data     = ft_appenddata([], data{:});
data     = removefields(data, {'elec' 'cfg'});

if mccaflag
  data = rt_mcca(data);
end

switch type
  case 'timelock'
      
    % Correct for the position of trials
    if poscorrect
        disp('Correct ERF for positions of trials.')
        % Calculate the average for each positions irrespective of phase
        x = data.trialinfo(:,1);
        [~, x] = sort(x);

        for k = 1:3
          % Check this approach, tlckpos(1) is not equal to position 1
          tmpcfg        = [];
          tmpcfg.trials = x(k:3:length(x));
          tmpcfg.preproc.baselinewindow = [-0.1 0];
          tmpcfg.preproc.demean         = 'yes';
          tlckpos(k) = ft_timelockanalysis(tmpcfg, data);
        end

        for k = 1:3
          tlckpos(k).time = tlckpos(1).time;
        end
        
        % Subtract average
        for k = 1:3
            k
            currtrials = x(k:3:length(x));
            for c = 1:length(currtrials)
                data.trial{currtrials(c)} = data.trial{currtrials(c)} - tlckpos(k).avg;
            end
        end
        
        
    end
    
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
      
      fname = fullfile(savepath, sprintf('%s_%s%s', subj, savename, suff));
      save(fname, 'tlck');
      clear data
    end
    
  case 'timelock123'
    x = data.trialinfo(:,1);
    [~, x] = sort(x);
    [~, x] = sort(x); % x now maps the first column's values onto 1:ntrl
    
    for k = 1:3
      tmpcfg        = [];
      tmpcfg.trials = find(mod(x-1,3)+1==k);
      tmpcfg.preproc.baselinewindow = [-0.1 0];
      tmpcfg.preproc.demean         = 'yes';
      tlck(k) = ft_timelockanalysis(tmpcfg, data);
    end
    
    for k = 1:3
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
      
      fname = fullfile(savepath, sprintf('%s_123%s%s', subj, savename, suff));
      save(fname, 'tlck');
      clear data
    end
    
  case 'tfr'
    if isempty(cfgfreq)
      cfgfreq = [];
      cfgfreq.method = 'mtmconvol';
      cfgfreq.foi    = 2:2:30;
      cfgfreq.t_ftimwin = ones(1,numel(cfgfreq.foi)).*0.5;
      cfgfreq.taper  = 'hanning';
      cfgfreq.toi    = (-150:15:450)./300;
      cfgfreq.pad    = 4;
    end
    
    conds = [1 2 3 5 6 7];
    for k = 1:numel(conds)
      cfgfreq.trials = find(data.trialinfo(:,2)==conds(k));
      if ~mccaflag
        freq(k)    = ft_combineplanar([], ft_freqanalysis(cfgfreq, data));
      else
        freq(k)    = ft_freqanalysis(cfgfreq, data);
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
      
      fname = fullfile(savepath, sprintf('%s_%s%s', subj, savename, suff));
      save(fname, 'freq');
      clear data
    end
    
    
  case 'rsa'
    % to-be-implemented
end

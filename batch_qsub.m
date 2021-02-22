% script to run a pipeline in a batch

savedir = '/project/3012026.13/jansch';

subjlist = {'sub-004' 'sub-005' 'sub-006' 'sub-008' 'sub-009' 'sub-010' ...
 'sub-011' 'sub-012' 'sub-014' 'sub-017' 'sub-018' 'sub-019' ...
 'sub-020' 'sub-022' 'sub-023' 'sub-024' 'sub-026' 'sub-027' ...
 'sub-028' 'sub-029' 'sub-031' 'sub-032' 'sub-033' 'sub-034' ...
 'sub-035' 'sub-037'};


cfgfreq = [];
cfgfreq.method = 'mtmconvol';
cfgfreq.foi    = 12.5:2.5:40;
cfgfreq.t_ftimwin = ones(1,numel(cfgfreq.foi)).*0.4;
cfgfreq.taper  = 'dpss';
cfgfreq.tapsmofrq = ones(1,numel(cfgfreq.foi)).*5;
cfgfreq.toi    = (-150:15:450)./300;
cfgfreq.pad    = 4;

for k = 1:numel(subjlist)
  % old syntax
  %qsubfeval('rt_sensorlevelanalysis', subjlist{k}, [], 'tfr', 1, savedir, 'memreq', 12*1024^3, 'timreq', 60*60, 'batchid', subjlist{k});

  qsubfeval('rt_sensorlevelanalysis', subjlist{k}, 'type', 'tfr', 'saveflag', savedir, 'savename', 'tfr_mid', 'cfgfreq', cfgfreq, 'memreq', 12*1024^3, 'timreq', 60*60, 'batchid', subjlist{k});

end

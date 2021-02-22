% script to run a pipeline in a batch

savedir = '/project/3012026.13/jansch';

subjlist = {'sub-004' 'sub-005' 'sub-006' 'sub-008' 'sub-009' 'sub-010' ...
 'sub-011' 'sub-012' 'sub-014' 'sub-017' 'sub-018' 'sub-019' ...
 'sub-020' 'sub-022' 'sub-023' 'sub-024' 'sub-026' 'sub-027' ...
 'sub-028' 'sub-029' 'sub-031' 'sub-032' 'sub-033' 'sub-034' ...
 'sub-035' 'sub-037'};

for k = 1:numel(subjlist)
  qsubfeval('rt_sensorlevelanalysis', subjlist{k}, [], 'timelock', 0, savedir, 'memreq', 12*1024^3, 'timreq', 60*60, 'batchid', subjlist{k});
end

function [atlasscaled] = rt_scaleatlas(cfgdata, atlas)

% This function will scale a provided atlas to the dimensions that are
% indicated by the cfg. The controlled scaling of the atlas guarantees that
% interpolation of individual parcels will be successful.

dim             = cfgdata.dim;

% Downsample atlas to get close to requested dimensions
atlasdim        = atlas.dim;
dimappr         = atlasdim./dim;
dimappr_high    = min(floor(dimappr));
dimappr_low     = min(ceil(dimappr));

cfg              = [];
cfg.downsample   = dimappr_high;
atlassample_high = ft_volumedownsample(cfg, atlas);
cfg.downsample   = dimappr_low;
atlassample_low  = ft_volumedownsample(cfg, atlas);


% Check the information loss of the downsampling

orgparc = unique(atlas.parcellation);
splparc = unique(atlassample.parcellation);
infloss = size(orgparc, 1) ~= size(splparc, 1);

if infloss
    
    % If parcel information is lost due to downsampling, figure out which
    % parcels were lost and how many 'voxels' were involved for each parcel
    
    lostparc = zeros(size(orgparc, 1)-size(splparc, 1),2);
    lostparcidx = 1;
    for numparc = 1:size(orgparc,1)
        if ~ismember(orgparc(numparc), splparc)
            lostparc(lostparcidx,1) = orgparc(numparc);
            lostparcidx = lostparcidx +1;
        end
    end
    
    % Count number of point for each lost parcel in original atlas
    for losts = 1:size(lostparc,1)
        lostparc(losts,2) = size(find(atlas.parcellation == lostparc(losts,1)),1);
    end
    
    disp("Lost parcels and #voxels are:")
    disp(lostparc)
    
else
    disp('No parcel information was lost due to downsampling.')
end
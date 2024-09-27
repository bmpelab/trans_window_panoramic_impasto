% apply 2d mask to 3d layers
% MASK : 2d mask
% LAYERS : 3d layers

% MASKED_LAYERS == 0, if MASK == true

%%
function [MASKED_LAYERS] = maskNdWith2d(MASK,LAYERS)

MASKED_LAYERS = LAYERS;

for i = 1 : size(LAYERS,3)
    layer = LAYERS(:,:,i);
    layer(MASK) = 0;
    MASKED_LAYERS(:,:,i) = layer;
end

end
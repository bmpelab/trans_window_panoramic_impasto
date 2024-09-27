function [SMASK] = smoothMaskGenerator(MASK1,MASK2)

H1 = size(MASK1,1);
W1 = size(MASK1,2);
H2 = size(MASK2,1);
W2 = size(MASK2,2);

dH = round((H1-H2)/2);
dW = round((W1-W2)/2);

se = strel('square',3);

mask2in1 = true(size(MASK1));
mask2in1(dH+1:dH+H2,dW+1:dW+W2) = MASK2;
mask_inner = imdilate(mask2in1,se);
mask_inner = mask_inner~=mask2in1;
mask_outside = imerode(mask2in1,se);
mask_outside = mask_outside~=mask2in1;
smask = logical(mask_outside+mask_inner);

mask_boundary = false(size(MASK1));
mask_boundary(dH+1:dH+H2,dW+1:dW+W2) = true;
mask_inner = imdilate(mask_boundary,se);
mask_inner = mask_inner~=mask_boundary;
mask_outside = imerode(mask_boundary,se);
mask_outside = mask_outside~=mask_boundary;
smask = logical(mask_outside+mask_inner);

SMASK = smask;

end
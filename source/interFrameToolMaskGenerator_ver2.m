% SFMAP : scene flow map (frame a->b)
% P3MAP : point 3d map (a)
% TPMASKA, TPMASKB : tool and point mask (a) and (b), == true means occluded
% PL : left camera projection matrix
% SR : sampling rate
% OFFSET : image offset

%%
function [IFTMASK] = interFrameToolMaskGenerator_ver2(SFMAP,P3MAP,TPMASKA,TPMASKB,PL,SR,OFFSET)

W = size(SFMAP,2);
H = size(SFMAP,1);
N = W*H;

flag_inter = true(N,1);

indices_a = 1 : N;
flag_tmaska = TPMASKA(:); % reshape(TMASKA,[],1);
points_a = reshape(P3MAP, [], 3);
sf = reshape(SFMAP, [], 3);

indices_b_in_a = indices_a(~flag_tmaska);
points_b = points_a(~flag_tmaska,:) + sf(~flag_tmaska,:);
% projection
pixels_b = PL*[points_b';ones(1,size(points_b,1))];
pixels_b = pixels_b./pixels_b(3,:);
pixels_b(1,:) = pixels_b(1,:)-OFFSET;
pixels_b = pixels_b/SR;
pixels_b = pixels_b(1:2,:)';
crs_b = round(pixels_b+1);
% boundary flag
flag_boundary = logical((crs_b(:,1)>=1).*(crs_b(:,1)<=W).*(crs_b(:,2)>=1).*(crs_b(:,2)<=H)); % within the boundary is true
indices_b_in_a = indices_b_in_a(flag_boundary);
%
indices_b = (crs_b(flag_boundary,1)-1)*H + crs_b(flag_boundary,2);
flag_tmaskb = TPMASKB(indices_b);
indices_b_in_a = indices_b_in_a(~flag_tmaskb);

flag_inter(indices_b_in_a) = false;

IFTMASK = reshape(flag_inter,H,W);

end
% SFMAP : h by w by 3 scene flow map
% DDMASK :h by w direct deformation mask from directDeformationMask
% P3MAP : h by w by 3 point 3d map
% IMAP : h by w index map
% N : number of vertices in hypermap

% I : m by n selection matrix, n is the size of points in the merged mesh,
    % m is the number of points with scene flow
% VD : vertices in desitination, VD = I*VS+SF
%%
function [I,VD,SF] = directDeformation_ver2(HYPERMAP,SFMAP,DDMASK,CCMASK)

H = size(SFMAP,1);
W = size(SFMAP,2);
HC = size(HYPERMAP.vertex_index_map,1);
WC = size(HYPERMAP.vertex_index_map,2);
dH = round((HC-H)/2);
dW = round((WC-W)/2);

IMAP = HYPERMAP.vertex_index_map(dH+1:dH+H,dW+1:dW+W);
IMAP(CCMASK(dH+1:dH+H,dW+1:dW+W)) = 0;
N = sum(IMAP(:)~=0);
IMAP(IMAP(:)~=0) = 1:N;
% IMAP = HYPERMAP.vertex_index_map;
% IMAP(CCMASK) = 0;
% N = sum(IMAP(:)~=0);
% IMAP(IMAP(:)~=0) = 1:N;
P3MAP = HYPERMAP.point_3d_map(dH+1:dH+H,dW+1:dW+W,:);

[I,VD,SF] = directDeformation(SFMAP,DDMASK,P3MAP,IMAP,N);

end
%%
function [I,VD,SF] = directDeformation(SFMAP,DDMASK,P3MAP,IMAP,N)

oflag = ~DDMASK(:); % reshape(~DDMASK,[],1); == true means valid scene flow

sf = reshape(SFMAP,[],3);
sf = sf(oflag,:);
SF = sf;

indices = IMAP(~DDMASK);
M = length(indices);

rows = (1:M)';
cols = indices(:);
values = ones(M,1);
I = sparse(rows,cols,values,M,N);

points = reshape(P3MAP,[],3);
points = points(oflag,:);
VD = SF + points;

end
% SFMAP : scene flow map (frame a->b)
% P3MAP : point 3d map (a)
% MASK = EDGEMASK+POINTMASK : EDGEMASK from edgeMaskGenerator or
    % surfaceReconstructionFromPoint3DMap, POINTMASK == true if point is (0,0,0)
% T : threshold

% MPSMASK : maximal principal strain mask, == true means mps higher than T or outlier

%%
function [MPSMASK,MPSMAP] = maxPrincipalStrainMaskGenerator(SFMAP,P3MAP,MASK,T)

W = size(SFMAP,2);
H = size(SFMAP,1);
N = W*H;

% scene flow vector : (u,v,w)
% 3d point : (x,y,z)
% displacement gradient : (du/dx du/dy du/dz dv/dx dv/dy dv/dz dw/dx dw/dy dw/dz)
mps_mask = true(H,W);
mps_map = zeros(H,W);
parfor i = 1 : N
    c = ceil(i/H);
    r = mod(i-1,H)+1;

    if MASK(r,c)
        continue;
    end
    
    ra = r-2;
    rb = r+2;
    if ra < 1
        ra = 1;
        rb = 4;
    elseif rb > H
        rb = H;
        ra = H-3;
    end

    ca = c-2;
    cb = c+2;
    if ca < 1
        ca = 1;
        cb = 4;
    elseif cb > W
        cb = W;
        ca = W-3;
    end

    neighbor_flags = reshape(MASK(ra:rb,ca:cb),[],1);
    neighbor_flags = ~neighbor_flags;
    if sum(neighbor_flags) < 13
        continue;
    end
    neighbor_flags((c-ca+1-1)*(rb-ra+1)+(r-ra+1)) = [];

    neighbor_points = reshape(P3MAP(ra:rb,ca:cb,:),[],3);
    central_point = neighbor_points((c-ca+1-1)*(rb-ra+1)+(r-ra+1),:); % P3MAP(r,c,:)
    neighbor_points((c-ca+1-1)*(rb-ra+1)+(r-ra+1),:) = []; % remove the point corresponding to (r,c)
    neighbor_points = neighbor_points(neighbor_flags,:);
    % coor is a num by 4 matrix
    % relative coordinates to central point
    % coor = |1 DeltaX1 DeltaY1 DeltaZ1|
    %        |1 DeltaX2 DeltaY2 DeltaZ2|
    %        |         ......          |
    coor = [ones(size(neighbor_points,1),1),neighbor_points - central_point];
    % coor*[displacement gradient] = disp
    % [displacement gradient] = dispGrad_smoothDisp is a 4 by 3 matrix
    % if rank(coor) == 4, which means coor has full column rank
    % its pseudo inverse matrix exists
    % note that coor_pi is the pseudo inverse matrix of coor
    % coor_pi is a 4 by num matrix
    % coor_pi*coor = I(4x4)
    % coor_pi*coor*dispGrad_smoothDisp = dispGrad_smoothDisp = coor_pi*disp
    if rank(coor) < 4 % coor do not have full column rank
        continue;
    end
    % calculate coor^(+)
    coor_pi = (coor'*coor)\coor'; % pi : pseudo inverse
    % disp is a num by 3 matrix
    % disp = |u1 v1 w1|
    %        |u2 v2 w2|
    %        | ...... |
    disp = reshape(SFMAP(ra:rb,ca:cb,:),[],3);
    disp((c-ca+1-1)*(rb-ra+1)+(r-ra+1),:) = [];
    disp = disp(neighbor_flags,:);
    % displacement gradient [4 by 3]
    % note: partial derivative := pd
    % dispGrad_smoothDisp = |     u0          v0          w0    |
    %                       |pd(u)/pd(x) pd(v)/pd(x) pd(w)/pd(x)|
    %                       |pd(u)/pd(y) pd(v)/pd(y) pd(w)/pd(y)|
    %                       |pd(u)/pd(z) pd(v)/pd(z) pd(w)/pd(z)|
    dispGrad_smoothDisp = coor_pi*disp;
    % strainTensor = [3 by 3] is a symmetric matrix
    %              = |epsilon_xx epsilon_xy epsilon_xz|
    %                |epsilon_yx epsilon_yy epsilon_yz|
    %                |epsilon_zx epsilon_zy epsilon_zz|
    strainTensor = zeros(3,3);
    strainTensor(1,1) = dispGrad_smoothDisp(2,1); % epsilon_xx = pd(u)/pd(x)
    strainTensor(2,2) = dispGrad_smoothDisp(3,2); % epsilon_yy = pd(v)/pd(y)
    strainTensor(3,3) = dispGrad_smoothDisp(4,3); % epsilon_zz = pd(w)/pd(z)
    strainTensor(1,2) = 0.5*(dispGrad_smoothDisp(3,1)+dispGrad_smoothDisp(2,2)); % epsilon_xy = 0.5*(pd(u)/pd(y)+pd(v)/pd(x))
    strainTensor(1,3) = 0.5*(dispGrad_smoothDisp(4,1)+dispGrad_smoothDisp(2,3)); % epsilon_xz = 0.5*(pd(u)/pd(z)+pd(w)/pd(x))
    strainTensor(2,3) = 0.5*(dispGrad_smoothDisp(4,2)+dispGrad_smoothDisp(3,3)); % epsilon_yz = 0.5*(pd(v)/pd(z)+pd(w)/pd(y))
    strainTensor(2,1) = strainTensor(1,2);
    strainTensor(3,1) = strainTensor(1,3);
    strainTensor(3,2) = strainTensor(2,3);
    % calculate the eigenvalues of strain tensor
    % using singular value decomposition (SVD)
    % principalStrain = sqrt(svd(strainTensor));
    % using eigen value decomposition (EVD) 
    principalStrain = eig(strainTensor);
    mps = max(principalStrain);
    mps_map(i) = mps;
    % the maximum principal strain
    if mps < T
        mps_mask(i) = false;
    end
end

MPSMASK = mps_mask;
MPSMAP = mps_map;

end
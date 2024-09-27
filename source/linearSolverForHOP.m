% linear solver for homeomorphism optimization problem

% CMASK : h by w constraint mask, == true is apply constraint, CMASK =
    % tool_mask.*~point_mask, means occluded and defined
% CMAP : h by w constraint map
% PL : left projection matrix
% SR : sampling rate
% OFFSET : image offset
% CF : constraint flag

% A*X=B
% X = A^-1*B

%%
function [X] = linearSolverForHOP(A,B,CMASK,CMAP,PL,SR,OFFSET,CF)

dA = decomposition(A);
X = dA\B;

if ~CF
    return;
end

W = size(CMASK,2);
H = size(CMASK,1);
N = size(X,1);

% project X to image plane at x
x = PL*[X';ones(1,N)];
x = x./x(3,:);
x(1,:) = x(1,:)-OFFSET;
x = x/SR;
x = x(1:2,:)';
crs = round(x+1);
% check whether x is occluded
% occluded is defined as: within the field of view, tool_mask = true
flag_fov = logical((crs(:,1)>=1).*(crs(:,1)<=W).*(crs(:,2)>=1).*(crs(:,2)<=H));
indices = (crs(:,1)-1)*H+crs(:,2);
flag_occluded = false(length(flag_fov),1);
flag_occluded(flag_fov) = CMASK(indices(flag_fov));
M = sum(flag_occluded);
% if non-occluded
if M == 0
    return;
end
% if occluded, nearest interpolate a constraint from CMAP
constraints = reshape(CMAP,[],3);
constraints = constraints(indices(flag_occluded),:);
% remove irrational constraint
flag_ir = logical(((constraints(:,3)-X(flag_occluded,3))<5));
flag_occluded(flag_occluded) = flag_ir;
M = sum(flag_occluded);
% if all irrational
if M == 0
    return;
end
%
constraints = constraints(flag_ir,:);
indices = find(flag_occluded);
rows = (1:M)';
cols = indices(:);
values = ones(M,1);
I = sparse(rows,cols,values,M,N);
%
X_archive = X;
X(:,3) = lsqlin(A,B(:,3),-I,-constraints(:,3));
directions = X_archive./sqrt(sum(X_archive.^2,2));
scales = (X(:,3) - X_archive(:,3))./directions(:,3);
X = X_archive + scales.*directions;

end
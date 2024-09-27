% MESH : triangulation representation of a mesh, mesh.Points, mesh.ConnectivityList
% PL : 3x4 left camera projection matrix
% W, H : width and height of the image
% SR : sample rate, this define the sub-pixel threshold

% FLAG : per-vertex flag, == true means visible to the camera; == false
    % means invisible to the camera, perhaps due to occlusion or out of the
    % field of view

%%
function [FLAG] = vertexVisibility(MESH,PL,W,H,SR)

num_point = size(MESH.Points,1);
depth_map = zeros(ceil(H/SR),ceil(W/SR));
index_map = zeros(size(depth_map));
pixels = PL*[MESH.Points';ones(1,num_point)];
pixels = pixels./pixels(3,:)/SR;
crs = round(pixels(1:2,:)'+1);
flag = true(num_point,1);
for i = 1 : num_point
    c = crs(i,1);
    r = crs(i,2);
    % negative z
    if MESH.Points(i,3) < 1e-6
        flag(i) = false;
        continue;
    end
    % out of the field of view
    if c < 1 || c > W || r < 1 || r > H
        flag(i) = false;
        continue;
    end
    if index_map(r,c) == 0
        index_map(r,c) = i;
        depth_map(r,c) = MESH.Points(i,3);
    else 
        % the current point is shallower than the stored point
        if depth_map(r,c) > MESH.Points(i,3)
            flag(index_map(r,c)) = false;
            index_map(r,c) = i;
            depth_map(r,c) = MESH.Points(i,3);
        % the current point is deeper than the stored point
        else
            flag(i) = false;
        end
    end
end
FLAG = flag;

end
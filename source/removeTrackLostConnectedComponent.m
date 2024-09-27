% MASK : logical(tool_mask_i+strain_mask) == true means lost of tracking
% for a connected component in hypermap.vertex_inedx_map, 
    % if the tracking of all of its elements gets lost, we remove this
    % connected components, otherwise the optimization method would fail
    % due to rank deficient in [im;alpha*dem]
%%
function [HYPERMAP2,FACES2,FOVMAP2] = removeTrackLostConnectedComponent(HYPERMAP,FACES,FOVMAP,MASK)

H1 = size(MASK,1);
W1 = size(MASK,2);
H2 = size(HYPERMAP.vertex_index_map,1);
W2 = size(HYPERMAP.vertex_index_map,2);
dH = round((H2-H1)/2);
dW = round((W2-W1)/2);

% generate the connected components of hypermap.vertex_inedx_map
ccmap = bwlabel(logical(HYPERMAP.vertex_index_map),4);
ccflag = false(max(ccmap(:)),1);
for r = 1 : H1
    for c = 1 : W1
        % lost of track
        if MASK(r,c)
            continue;
        end
        % in track
        % check the cc label
        cclabel = ccmap(r+dH,c+dW);
        % the background
        if cclabel == 0
            continue;
        end
        %
        ccflag(cclabel) = true;
    end
end
if sum(ccflag==true)==length(ccflag)
    HYPERMAP2 = HYPERMAP;
    FACES2 = FACES;
    FOVMAP2 = FOVMAP;
    return;
end
% remove info in hypermap and faces
face_flag = true(size(FACES,1),1);
vertex_flag = true(HYPERMAP.num_of_vertex,1);
hypermap2 = HYPERMAP;
fovmap2 = FOVMAP;
for r = 1 : H2
    for c = 1 : W2
        vertex_id = HYPERMAP.vertex_index_map(r,c);
        if vertex_id == 0
            continue;
        end
        cclabel = ccmap(r,c);
        if ccflag(cclabel)
            continue;
        end
        %
        vertex_flag(vertex_id) = false;
        %
        face_ids = squeeze(FOVMAP(r,c,:));
        face_ids = face_ids(face_ids~=0);
        face_flag(face_ids) = false;
        %
        hypermap2.vertex_index_map(r,c) = 0;
        hypermap2.point_3d_map(r,c,:) = 0;
        hypermap2.derivative_map(r,c,:) = 0;
        hypermap2.derivative_direction_map(r,c,:) = -1;
        hypermap2.color_map(r,c,:) = uint8(0);
        hypermap2.backward_scene_flow_map(r,c,:) = 0;
        fovmap2(r,c,:) = 0;
    end
end
% update index map
num_valid_vertex = sum(vertex_flag);
hypermap2.num_of_vertex = num_valid_vertex;
vertex_indices = reshape(hypermap2.vertex_index_map,[],1);
vertex_indices(vertex_indices~=0) = 1 : num_valid_vertex;
hypermap2.vertex_index_map = reshape(vertex_indices, 2*H1, 2*W1);
% update the indices in faces
vertex_index_update = zeros(HYPERMAP.num_of_vertex,1);
vertex_index_update(vertex_flag) = 1 : num_valid_vertex;
faces2 = FACES(face_flag,:);
faces2(:) = vertex_index_update(faces2(:));
FACES2 = faces2;
% update the indices in fov_index_map
face_index_update = zeros(length(face_flag),1);
face_index_update(face_flag) = 1 : sum(face_flag);
fov_index_array = fovmap2(:);
fov_index_array(fov_index_array~=0) = face_index_update(fov_index_array(fov_index_array~=0));
fovmap2(:) = fov_index_array;
FOVMAP2 = fovmap2;

HYPERMAP2 = hypermap2;

end
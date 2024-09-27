%% ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ DEFORMATION RECOVERY ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
module_path = './source';
addpath(module_path);

texture_flag = true;
save_flag = true;
edge_filter_threshold = 10;
edge_filter_threshold2 = 1;
sr = 3; % sampling rate
alpha = 1.07;
image_offset = 64;
camera_fix_flag = true;
data_folder = './surgem_ex_vivo/g2';
left_image_folder = [data_folder '/rectified_left'];
tool_mask_folder = [data_folder '/mask'];
point_3d_map_folder = [data_folder '/point_3d_map'];
scene_flow_map_folder = [data_folder '/scene_flow'];
constraint_map_folder = [data_folder '/constraint_map'];
load([data_folder '/rectifiedCamera.mat']);
% point_3d_map
point_3d_map_dir = dir([point_3d_map_folder '/*.mat']);
point_3d_map_names = {point_3d_map_dir.name};
[point_3d_map_names,~] = sortNat(point_3d_map_names);
% left image
left_image_dir = dir([left_image_folder '/*.png']);
left_image_names = {left_image_dir.name};
[left_image_names,~] = sortNat(left_image_names);
% left mask
tool_mask_dir = dir([tool_mask_folder '/*.png']);
tool_mask_names = {tool_mask_dir.name};
[tool_mask_names,~] = sortNat(tool_mask_names);
% scene flow file (.mat) list
scene_flow_map_dir = dir([scene_flow_map_folder '/*.mat']);
scene_flow_map_names = {scene_flow_map_dir.name};
[scene_flow_map_names,~] = sortNat(scene_flow_map_names);
% constraint map file (.mat) list
constraint_map_dir = dir([constraint_map_folder '/*.mat']);
constraint_map_names = {constraint_map_dir.name};
[constraint_map_names,~] = sortNat(constraint_map_names);
% output mesh folder
texture_mesh_folder = [data_folder '/mesh'];
if ~exist(texture_mesh_folder,'dir')
    mkdir(texture_mesh_folder);
end

left_image = imread([left_image_folder '/' left_image_names{1}]);
left_image = left_image(:,1+image_offset:end,:);
left_image = left_image(1:sr:end,1:sr:end,:);
point_3d_map = load([point_3d_map_folder '/' point_3d_map_names{1}]).point_3d_map;
point_3d_map = point_3d_map(:,1+image_offset:end,:);
point_3d_map = point_3d_map(1:sr:end,1:sr:end,:);
point_mask = logical((point_3d_map(:,:,1)==0).*(point_3d_map(:,:,2)==0).*(point_3d_map(:,:,3)==0));
tool_mask = imread([tool_mask_folder '/' tool_mask_names{1}]);
tool_mask = tool_mask(:,1+image_offset:end,:);
tool_mask = tool_mask(1:sr:end,1:sr:end,:);
tool_mask = logical((tool_mask(:,:,1)==0).*(tool_mask(:,:,2)==255).*(tool_mask(:,:,3)==0));
% tool_mask = logical(tool_mask);
height = size(left_image,1);
width = size(left_image,2);
heightc = round(1.2*height); % canonical
widthc = round(1.2*width);
dh = round((heightc-height)/2);
dw = round((widthc-width)/2);

%
[edge_mask,vertex_index_map,fov_index_map,vertices,hedges,vedges,dedges,faces,colors] = surfaceReconstructionFromPoint3DMap_ver6...
    (point_3d_map,logical(tool_mask+point_mask),left_image,1,false,true,edge_filter_threshold,height*width/100);
derivative_map = derivativeMapGenerator_ver2(point_3d_map,~logical(vertex_index_map));

%
hypermap.color_map = uint8(zeros(heightc,widthc,3));
hypermap.color_map(dh+1:dh+height,dw+1:dw+width,:) = maskNdWith2d(~logical(vertex_index_map),left_image);
hypermap.backward_scene_flow_map = zeros(heightc,widthc,3);
hypermap.forward_scene_flow_map = zeros(heightc,widthc,3);
hypermap.backward_optical_flow_map = zeros(heightc,widthc,2);
hypermap.forward_optical_flow_map = zeros(heightc,widthc,2);
hypermap.point_3d_map = zeros(heightc,widthc,3);
hypermap.point_3d_map(dh+1:dh+height,dw+1:dw+width,:) = maskNdWith2d(~logical(vertex_index_map),point_3d_map);
hypermap.vertex_index_map = zeros(heightc,widthc);
hypermap.vertex_index_map(dh+1:dh+height,dw+1:dw+width) = vertex_index_map;
hypermap.num_of_vertex = max(vertex_index_map(:));
edge_mask_hyper = true(heightc,widthc);
edge_mask_hyper(dh+1:dh+height,dw+1:dw+width) = edge_mask;
[hypermap.derivative_map,hypermap.derivative_direction_map] = derivativeMapGenerator_ver3(hypermap.point_3d_map,~logical(hypermap.vertex_index_map),edge_mask_hyper);
fov_index_map_hyper = zeros(heightc,widthc,6);
fov_index_map_hyper(dh+1:dh+height,dw+1:dw+width,:) = maskNdWith2d(~logical(vertex_index_map),fov_index_map);
faces_hyper = faces;
bindex_map = zeros(size(hypermap.vertex_index_map));
%
point_3d_map_b = point_3d_map;
derivative_map_b = derivative_map;
point_mask_b= point_mask;
tool_mask_b = tool_mask;
edge_mask_b = edge_mask;
vertex_index_map_b = vertex_index_map;
%
for i = 1 : length(scene_flow_map_names)
    disp(i);
    % left image (color map)
    left_image_b = imread([left_image_folder '/' left_image_names{i+1}]);
    left_image_b = left_image_b(:,1+image_offset:end,:);
    left_image_b = left_image_b(1:sr:end,1:sr:end,:);
    % load scene flow map
    scene_flow_map = load([scene_flow_map_folder '/' scene_flow_map_names{i}]).scene_flow_map; % a -> b
    scene_flow_map = scene_flow_map(:,1+image_offset:end,:);
    scene_flow_map = scene_flow_map(1:sr:end,1:sr:end,:);
    scene_flow_mask = logical((scene_flow_map(:,:,1)==0).*(scene_flow_map(:,:,2)==0).*(scene_flow_map(:,:,3)==0));
    % tool mask
    tool_mask_a = tool_mask_b;
    tool_mask_b = imread([tool_mask_folder '/' tool_mask_names{i+1}]);
    tool_mask_b = tool_mask_b(:,1+image_offset:end,:);
    tool_mask_b = tool_mask_b(1:sr:end,1:sr:end,:);
    tool_mask_b = logical((tool_mask_b(:,:,1)==0).*(tool_mask_b(:,:,2)==255).*(tool_mask_b(:,:,3)==0));
    % point 3d map
    point_3d_map_a = point_3d_map_b;
    point_3d_map_b = load([point_3d_map_folder '/' point_3d_map_names{i+1}]).point_3d_map;
    point_3d_map_b = point_3d_map_b(:,1+image_offset:end,:);
    point_3d_map_b = point_3d_map_b(1:sr:end,1:sr:end,:);
    point_mask_a = point_mask_b;
    point_mask_b = logical((point_3d_map_b(:,:,1)==0).*(point_3d_map_b(:,:,2)==0).*(point_3d_map_b(:,:,3)==0));
    %
    vertex_index_map_a = vertex_index_map_b;
    edge_mask_a = edge_mask_b;
    [edge_mask_b,vertex_index_map_b,~,~,~,~,~,~,~]...
        = surfaceReconstructionFromPoint3DMap_ver6(point_3d_map_b,logical(tool_mask_b+point_mask_b),[],1,false,true,edge_filter_threshold,height*width/100);
    % inter-frame tool mask
    [tool_mask_i] = interFrameToolMaskGenerator_ver2(scene_flow_map,point_3d_map_a,logical(~logical(vertex_index_map_a)+scene_flow_mask),...
        ~logical(vertex_index_map_b),PL,sr,image_offset);
    % strain-based filtering
    % undefined, discoutinuous and occluded area = ~logical(vertex_index_map_a)
    % undefined and discoutinuous area = ~logical(vertex_index_map_a).*~tool_mask_a
    [strain_mask,strain_map] = maxPrincipalStrainMaskGenerator(scene_flow_map,point_3d_map_a,...
        logical(~logical(vertex_index_map_a).*~tool_mask_a+scene_flow_mask),1);
    %
    derivative_map_a = derivative_map_b;
    [derivative_map_b,~] = derivativeMapGenerator_ver3(point_3d_map_b,~logical(vertex_index_map_b),edge_mask_b);
    % check the existance of valid scene flow
    [dd_mask] = directDeformationMask(scene_flow_map,logical(tool_mask_i+strain_mask),hypermap.vertex_index_map);
    %
    [cc_map,cc_mask,cc_group_id,cc_centers,iso_flags] = connectedComponentGrouping(dd_mask,hypermap.vertex_index_map);
    % per-vertex scene flow assignment
    [im,pt,sft] = directDeformation_ver2(hypermap,scene_flow_map,dd_mask,cc_mask); % im : selection matrix; pt : direct update of point
    %
    if ~camera_fix_flag
        [~,rotation,translation,~] = absoluteOrientationQuaternion((pt-sft)', pt');
    else
        rotation = eye(3);
        translation = zeros(3,1);
    end
    pose = [[rotation,translation];[0 0 0 1]];
    disp(pose);
    if ~camera_fix_flag
        % rotate the derivetives
        derivatives = reshape(hypermap.derivative_map,[],6);
        derivatives(:,1:3) = (rotation*derivatives(:,1:3)')';
        derivatives(:,4:6) = (rotation*derivatives(:,4:6)')';
        hypermap.derivative_map = reshape(derivatives,heightc,widthc,[]);
    end
    %
    [dem,dce] = dM_ver5(hypermap,cc_mask,height,width);
    % homeomorphism optimization
    % optimization variable: M_pe: U_ps -> S_pe
    constraint_map = load([constraint_map_folder '/' constraint_map_names{i+1}]).constraint_map;
    constraint_map = constraint_map(:,1+image_offset:end,:);
    constraint_map = constraint_map(1:sr:end,1:sr:end,:);
    % if there exist a connected component without any direct deformation information,
        % the optimization would fail due to rank deficient in [im;alpha*dem]
    constraint_mask = tool_mask_b;
    [pe] = linearSolverForHOP([im;alpha*dem],[pt;alpha*dce],constraint_mask,constraint_map,PL,sr,image_offset,true);
    % save the forward scene flow map into hypermap_a
    [hypermap,pe,forward_scene_flow_hyper,forward_optical_flow_hyper]...
        = saveForwardFlow(hypermap,pe,cc_map,cc_mask,cc_group_id,cc_centers,iso_flags,alpha,pose,PL,sr,image_offset);
    % save mesh_a
    if save_flag
        if i == 1
            hypermap_output = hypermap;
        end
    end
    if texture_flag && save_flag
        [mesh] = surfaceReconstructionFromHypermap(hypermap_output,texture_flag);
        writeSurfaceMesh(mesh,[texture_mesh_folder '\' erase(point_3d_map_names{i},'mat') 'ply']);
    end
    % local parameterization of S_pe by projecting onto the canonical image plane (2h by 2w)
    % only keep the invisible part as S_po
    % use the index map for reducing the searching space for acceleration
    % U_po -> Spo : point 3d map (2h by 2w by 3)
    % U_po -> Sao : anchor 3d map (last visible point 3d map) (2h by 2w by 3)
    index_flag_hyper = reshape(hypermap.vertex_index_map~=0,[],1);
    colors_hyper = reshape(hypermap.color_map,[],3);
    colors_hyper = colors_hyper(index_flag_hyper,:);
    derivatives_hyper = reshape(hypermap.derivative_map,[],6);
    derivatives_hyper = derivatives_hyper(index_flag_hyper,:);
    bindices = bindex_map(:);
    bindices = bindices(index_flag_hyper);
    % fill hole if necessary
    % para_mask == true, means no information
    [point_3d_map_para,derivative_map_para,color_map_para,backward_scene_flow_map_para,backward_optical_flow_map_para,mask_para,bindex_map] = ...
        parameterization_ver7(height,width,heightc,widthc,pe,faces_hyper,derivatives_hyper,colors_hyper,...
        -forward_scene_flow_hyper,-forward_optical_flow_hyper,bindices,PL,sr,image_offset);
    % merge
    [point_3d_map_hyper,derivative_map_hyper,color_map_hyper,point_mask_hyper,bindex_map]...
        = mergeMap_ver2(point_3d_map_para,point_3d_map_b,...
        derivative_map_para,derivative_map_b,...
        color_map_para,left_image_b,...
        mask_para,~logical(vertex_index_map_b),tool_mask_b,bindex_map,false,[]);
  
    if save_flag
        hypermap_output = hypermap;
        [edge_mask_output,vertex_index_map_output,~,~,~,~,~,~,~]...
            = surfaceReconstructionFromPoint3DMap_ver6(point_3d_map_hyper,point_mask_hyper,[],1,false,true,edge_filter_threshold2*edge_filter_threshold,height*width/100);
        hypermap_output.point_3d_map = maskNdWith2d(~logical(vertex_index_map_output),point_3d_map_hyper);
        hypermap_output.color_map = maskNdWith2d(~logical(vertex_index_map_output),color_map_hyper);
        hypermap_output.backward_scene_flow_map = maskNdWith2d(~logical(vertex_index_map_output),backward_scene_flow_map_para);
        hypermap_output.backward_optical_flow_map = maskNdWith2d(~logical(vertex_index_map_output),backward_optical_flow_map_para);
        hypermap_output.derivative_map = maskNdWith2d(~logical(vertex_index_map_output),derivative_map_hyper);
        hypermap_output.derivative_direction_map = derivativeDirectionMapGenerator(~logical(vertex_index_map_output),edge_mask_output);
        hypermap_output.vertex_index_map = vertex_index_map_output;
        hypermap_output.num_of_vertex = max(vertex_index_map_output(:));
    end

    % apply continuity filtering or not
    [edge_mask_hyper,vertex_index_map_hyper,fov_index_map_hyper,~,~,~,~,faces_hyper,~]...
        = surfaceReconstructionFromPoint3DMap_ver6(point_3d_map_hyper,point_mask_hyper,[],1,false,true,edge_filter_threshold2*edge_filter_threshold,height*width/100);
    %
    hypermap.point_3d_map = maskNdWith2d(~logical(vertex_index_map_hyper),point_3d_map_hyper);
    hypermap.color_map = maskNdWith2d(~logical(vertex_index_map_hyper),color_map_hyper);
    hypermap.backward_scene_flow_map = maskNdWith2d(~logical(vertex_index_map_hyper),backward_scene_flow_map_para);
    hypermap.backward_optical_flow_map = maskNdWith2d(~logical(vertex_index_map_hyper),backward_optical_flow_map_para);
    hypermap.derivative_map = maskNdWith2d(~logical(vertex_index_map_hyper),derivative_map_hyper);
    hypermap.derivative_direction_map = derivativeDirectionMapGenerator(~logical(vertex_index_map_hyper),edge_mask_hyper);
    hypermap.vertex_index_map = vertex_index_map_hyper;
    hypermap.num_of_vertex = max(vertex_index_map_hyper(:));

    if i == length(scene_flow_map_names)
        % save hypermap without forward scene flow map
        % save mesh_b
        if texture_flag && save_flag
            [mesh] = surfaceReconstructionFromHypermap(hypermap_output,texture_flag);
            writeSurfaceMesh(mesh,[texture_mesh_folder '\' erase(point_3d_map_names{i+1},'mat') 'ply']);
        end
    end
end

%%
function parsave(path,var_name,var)

rename_var(var_name,var)
save(path,var_name);

end
%
function rename_var(var_name,var)

assignin('caller',var_name,var);

end
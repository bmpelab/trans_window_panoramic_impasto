% PE : [pe] = linearSolverForHOP([im;alpha*dem],[pt;alpha*dce],tool_mask_b,point_3d_map_b,PL,sr);

%%
function [HYPERMAP2,SFS] = saveForwardSceneFlow(HYPERMAP,PE)

index_flag = reshape(HYPERMAP.vertex_index_map~=0,[],1);
ps = reshape(HYPERMAP.point_3d_map,[],3); % source point
forward_scene_flow = reshape(HYPERMAP.forward_scene_flow_map,[],3);
SFS = PE-ps(index_flag,:);
forward_scene_flow(index_flag,:) = SFS;
HYPERMAP.forward_scene_flow_map = reshape(forward_scene_flow,size(HYPERMAP.forward_scene_flow_map,1),size(HYPERMAP.forward_scene_flow_map,2),[]);

HYPERMAP2 = HYPERMAP;

end
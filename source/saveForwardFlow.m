% PE : [pe] = linearSolverForHOP([im;alpha*dem],[pt;alpha*dce],tool_mask_b,point_3d_map_b,PL,sr);
% CCMASK : == true means no pe

%%
function [HYPERMAP2,PE2,SFS,OFS] = saveForwardFlow(HYPERMAP,PE,CC_MAP,CC_MASK,CC_GROUP_IDS,CC_CENTERS,ISO_FLAGS,ALPHA,POSE,PL,SR,OFFSET)

HC = size(CC_MAP,1);
WC = size(CC_MAP,2);

index_flag = reshape(~CC_MASK,[],1);
ps = reshape(HYPERMAP.point_3d_map,[],3); % source point
ccids = reshape(CC_MAP,[],1);
forward_scene_flow = reshape(HYPERMAP.forward_scene_flow_map,[],3);
forward_scene_flow(index_flag,:) = PE-ps(index_flag,:);
forward_scene_flow_map = reshape(forward_scene_flow,HC,WC,[]);
for i = 1 : length(CC_GROUP_IDS)
    if CC_GROUP_IDS(i) == i
        continue;
    end
    if ISO_FLAGS(i)
        % isolated connected component
        % the scene flow of this group (i) is equal to the scene flow of the
            % closest element of the linked group
%         [rs,cs] = find(CC_MAP==CC_GROUP_IDS(i));
%         cc_center = CC_CENTERS(i,:); % (c,r)
%         distances = sum(([cs,rs]-cc_center).^2,2);
%         [~,min_i] = min(distances);
%         forward_scene_flow(ccids==i,:) = repmat(squeeze(forward_scene_flow_map(rs(min_i),cs(min_i),:))',sum(ccids==i),1);

%         forward_scene_flow(ccids==i,:) = repmat(mean(forward_scene_flow(ccids==CC_GROUP_IDS(i),:)),sum(ccids==i),1);

          ps_ = ps(ccids==i,:);
          pe_ = POSE*[ps_';ones(1,size(ps_,1))];
          pe_ = pe_(1:3,:)';
          forward_scene_flow(ccids==i,:) = pe_ - ps_;
    else
        se = strel("diamond",1);
        mask1 = (CC_MAP==i);
        mask3 = logical(~CC_MASK.*imdilate(mask1,se));
        index_map = mask1+mask3;
        n = sum(index_map(:)~=0);
        index_map(index_map(:)~=0) = 1 : n;
        index_mask1 = index_map(mask1);
        n1 = length(index_mask1);

        bool_map = false(3,3);
        bool_map(2,2) = true;
        bool_map(2,3) = true;
        se2 = strel(bool_map);
        bmask_lr = logical((logical(mask1+mask3)~=imerode(logical(mask1+mask3),se2)));
        bool_map = false(3,3);
        bool_map(2,2) = true;
        bool_map(3,2) = true;
        se2 = strel(bool_map);
        bmask_tb = logical((logical(mask1+mask3)~=imerode(logical(mask1+mask3),se2)));
        [dem,dce] = dM_ver7(HYPERMAP,~logical(mask1+mask3),bmask_lr,bmask_tb,POSE(1:3,1:3),1);
        [im,pt,sft] = directDeformation_ver2(HYPERMAP,forward_scene_flow_map,~mask3,~logical(mask1+mask3));

%         dce0 = zeros(size(dce));
%         A = [im;ALPHA*dem];
%         B = [sft;ALPHA*dce0];
%         dA = decomposition(A);
%         sfe = dA\B;
%         forward_scene_flow(ccids==i,:) = sfe(index_mask1,:);

        A = [im;ALPHA*dem];
        B = [pt;ALPHA*dce];
        dA = decomposition(A);
        pe = dA\B;
        forward_scene_flow(ccids==i,:) = pe(index_mask1,:) - ps(ccids==i,:);

%         im2 = sparse((1:n1)',index_mask1,ones(n1,1),n1,n);
%         pose = POSE;
%         pose(1:3,4) = mean(sft)';
%         pt2 = pose*[ps(ccids==i,:)';ones(1,n1)];
%         pt2 = pt2(1:3,:)';
%         A = [im;im2];
%         B = [pt;pt2];
%         dA = decomposition(A);
%         pe = dA\B;
%         forward_scene_flow(ccids==i,:) = pe(index_mask1,:) - ps(ccids==i,:);
    end
end

%
index_flag2 = reshape(logical(HYPERMAP.vertex_index_map),[],1);
pe2 = ps(index_flag2,:) + forward_scene_flow(index_flag2,:);
PE2 = pe2;

SFS = forward_scene_flow(index_flag2,:);

HYPERMAP.forward_scene_flow_map = reshape(forward_scene_flow,HC,WC,[]);

pixels_s = PL*[ps(index_flag2,:)';ones(1,HYPERMAP.num_of_vertex)];
pixels_s = pixels_s./pixels_s(3,:);
pixels_s = pixels_s(1:2,:)';
pixels_s(:,1) = pixels_s(:,1) - OFFSET;
pixels_s = pixels_s/SR;
pixels_e = PL*[pe2';ones(1,HYPERMAP.num_of_vertex)];
pixels_e = pixels_e./pixels_e(3,:);
pixels_e = pixels_e(1:2,:)';
pixels_e(:,1) = pixels_e(:,1) - OFFSET;
pixels_e = pixels_e/SR;
forward_optical_flow = reshape(HYPERMAP.forward_optical_flow_map,[],2);
forward_optical_flow(index_flag2,:) = pixels_e - pixels_s;
OFS = forward_optical_flow(index_flag2,:);
HYPERMAP.forward_optical_flow_map = reshape(forward_optical_flow,HC,WC,[]);

HYPERMAP2 = HYPERMAP;

end
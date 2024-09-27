% P3MAP : point 3d map

%%
function [EDGE_MASK] = edgeMaskGenerator(P3MAP)

height = size(P3MAP,1);
width = size(P3MAP,2);

edge_8_map = -ones(height,width,8);
edge_h_map = sqrt(sum((P3MAP(:,1:end-1,:)-P3MAP(:,2:end,:)).^2,3)); % horizontal
edge_8_map(:,1:end-1,1) = edge_h_map; % to the right
edge_8_map(:,2:end,2) = edge_h_map; % to the left
edge_v_map = sqrt(sum((P3MAP(1:end-1,:,:)-P3MAP(2:end,:,:)).^2,3)); % vertical
edge_8_map(1:end-1,:,3) = edge_v_map; % to the bottom
edge_8_map(2:end,:,4) = edge_v_map; % to the top
edge_d1_map = sqrt(sum((P3MAP(1:end-1,1:end-1,:)-P3MAP(2:end,2:end,:)).^2,3)); % diagonal \
edge_8_map(1:end-1,1:end-1,5) = edge_d1_map; % to the bottom-right
edge_8_map(2:end,2:end,6) = edge_d1_map; % to the top-left
edge_d2_map = sqrt(sum((P3MAP(2:end,1:end-1,:)-P3MAP(1:end-1,2:end,:)).^2,3)); % diagonal /
edge_8_map(2:end,1:end-1,7) = edge_d2_map; % to the top-right
edge_8_map(1:end-1,2:end,8) = edge_d2_map; % to the bottom-left
edge_m_map = sum(edge_8_map,3)./sum(edge_8_map~=-1,3); % mean
edge_m_map(edge_m_map<1e-6) = 1e-6;

edge_d8_map = -ones(height,width,8);
edge_dh_map = edge_m_map(:,1:end-1)./edge_m_map(:,2:end); % differential horizontal
edge_d8_map(:,1:end-1,1) = edge_dh_map; % to the right
edge_d8_map(:,2:end,2) = edge_dh_map; % to the left
edge_dv_map = edge_m_map(1:end-1,:)./edge_m_map(2:end,:); % differential vertical
edge_d8_map(1:end-1,:,3) = edge_dv_map; % to the bottom
edge_d8_map(2:end,:,4) = edge_dv_map; % to the top
edge_dd1_map = edge_m_map(1:end-1,1:end-1,:)./edge_m_map(2:end,2:end,:); % diagonal \
edge_d8_map(1:end-1,1:end-1,5) = edge_dd1_map; % to the bottom-right
edge_d8_map(2:end,2:end,6) = edge_dd1_map; % to the top-left
edge_dd2_map = edge_m_map(2:end,1:end-1,:)./edge_m_map(1:end-1,2:end,:); % diagonal /
edge_d8_map(2:end,1:end-1,7) = edge_dd2_map; % to the top-right
edge_d8_map(1:end-1,2:end,8) = edge_dd2_map; % to the bottom-left
edge_d8_map(logical((edge_d8_map<1).*(edge_d8_map>0))) = 1./edge_d8_map(logical((edge_d8_map<1).*(edge_d8_map>0)));
edge_dm_map = sum(edge_d8_map,3)./sum(edge_d8_map~=-1,3); % mean

edge_mask = edge_dm_map<1.5;
edge_mask = bwareafilt(edge_mask,[round(height*width/100) inf]);
edge_mask = ~edge_mask;
EDGE_MASK = edge_mask;

end
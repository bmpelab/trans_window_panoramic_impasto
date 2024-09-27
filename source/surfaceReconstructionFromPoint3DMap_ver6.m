% P3MAP : point 3d map
% MASK : if MASK(r,c) == true, ignore the vertex of the pixel (c-1,r-1)
% CMAP : color map (RGB image)
% SR : sample rate
% CF : color flag
% EF : edge-based filter flag
% ET : edge-based filter threshold
% ET2 : edge-based filter threshold

% V_INDEX_MAP : vertex index map
% FOV_INDEX_MAP : faces of vertex index map

%%
function [EDGE_MASK,V_INDEX_MAP,FOV_INDEX_MAP,VERTICES,HEDGES,VEDGES,DEDGES,FACES,COLORS]...
    = surfaceReconstructionFromPoint3DMap_ver6(P3MAP,MASK,CMAP,SR,CF,EF,ET,ET2)

height = size(P3MAP,1);
width = size(P3MAP,2);

p3map = P3MAP(1:SR:height,1:SR:width,:);
mask = MASK(1:SR:height,1:SR:width);
if CF
    cmap = CMAP(1:SR:height,1:SR:width,:);
else
    cmap = [];
end

[EDGE_MASK,V_INDEX_MAP,FOV_INDEX_MAP,VERTICES,HEDGES,VEDGES,DEDGES,FACES,COLORS] = surfaceReconstructionFromPoint3DMap(p3map,mask,cmap,CF,EF,ET,ET2);

end
%%
function [EDGE_MASK,V_INDEX_MAP,FOV_INDEX_MAP,VERTICES,HEDGES,VEDGES,DEDGES,FACES,COLORS]...
    = surfaceReconstructionFromPoint3DMap(P3MAP,MASK,CMAP,CF,EF,ET,ET2)

height = size(P3MAP,1);
width = size(P3MAP,2);

p1 = 10;
p2 = 90;
scale = ET;
outlier_8_map = false(height,width,8); % == false is outlier
edge_h_map = sqrt(sum((P3MAP(:,1:end-1,:)-P3MAP(:,2:end,:)).^2,3)); % horizontal
t1 = prctile(edge_h_map(edge_h_map(:)~=0),p1);
t2 = prctile(edge_h_map(edge_h_map(:)~=0),p2);
outlier_8_map(:,1:end-1,1) = logical(~MASK(:,1:end-1).*((edge_h_map<t2+scale*(t2-t1))+MASK(:,2:end))); % to the right
outlier_8_map(:,2:end,2) = logical(~MASK(:,2:end).*((edge_h_map<t2+scale*(t2-t1))+MASK(:,1:end-1))); % to the left
edge_v_map = sqrt(sum((P3MAP(1:end-1,:,:)-P3MAP(2:end,:,:)).^2,3)); % vertical
t1 = prctile(edge_v_map(edge_v_map(:)~=0),p1);
t2 = prctile(edge_v_map(edge_v_map(:)~=0),p2);
outlier_8_map(1:end-1,:,3) = logical(~MASK(1:end-1,:).*((edge_v_map<t2+scale*(t2-t1))+MASK(2:end,:))); % to the bottom
outlier_8_map(2:end,:,4) = logical(~MASK(2:end,:).*((edge_v_map<t2+scale*(t2-t1))+MASK(1:end-1,:))); % to the top
edge_d1_map = sqrt(sum((P3MAP(1:end-1,1:end-1,:)-P3MAP(2:end,2:end,:)).^2,3)); % diagonal \
t1 = prctile(edge_d1_map(edge_d1_map(:)~=0),p1);
t2 = prctile(edge_d1_map(edge_d1_map(:)~=0),p2);
outlier_8_map(1:end-1,1:end-1,5) = logical(~MASK(1:end-1,1:end-1).*((edge_d1_map<t2+scale*(t2-t1))+MASK(2:end,2:end))); % to the bottom-right
outlier_8_map(2:end,2:end,6) = logical(~MASK(2:end,2:end).*((edge_d1_map<t2+scale*(t2-t1))+MASK(1:end-1,1:end-1))); % to the top-left
edge_d2_map = sqrt(sum((P3MAP(2:end,1:end-1,:)-P3MAP(1:end-1,2:end,:)).^2,3)); % diagonal /
t1 = prctile(edge_d2_map(edge_d2_map(:)~=0),p1);
t2 = prctile(edge_d2_map(edge_d2_map(:)~=0),p2);
outlier_8_map(2:end,1:end-1,7) = logical(~MASK(2:end,1:end-1).*((edge_d2_map<t2+scale*(t2-t1))+MASK(1:end-1,2:end))); % to the top-right
outlier_8_map(1:end-1,2:end,8) = logical(~MASK(1:end-1,2:end).*((edge_d2_map<t2+scale*(t2-t1))+MASK(2:end,1:end-1))); % to the bottom-left

edge_mask = logical(prod(outlier_8_map,3));
edge_mask = bwareafilt(edge_mask,[ET2 inf]);
edge_mask = ~edge_mask;
EDGE_MASK = edge_mask;

if EF
    mask = logical(edge_mask+MASK);
else
    mask = MASK;
end
se = strel('square',3);
mask = imdilate(mask,se);
mask = imerode(mask,se);

vertex_flag = true(width*height,1);
face_flag = true(2*(width-1)*(height-1),1);
edge_flag = false((width-1)*height+width*(height-1)+(width-1)*(height-1),1); % horizontal + vertical + diagonal /
fov_index_map = zeros(height,width,6);
for c = 1 : width
    for r = 1 : height
        vertex_index = (c-1)*height+r;
        [face_indices,faces,edge_indices] = findNeighborFacesEdges(r,c,height,width);
        % check mask
        if mask(r,c)
            vertex_flag(vertex_index) = false;
            face_flag(face_indices) = false;
            continue;
        end

        % check mask of the neighbor
        neighbor_face_flag = true(length(face_indices),1);
        for i = 1 : length(neighbor_face_flag)
            for j = 2 :3
                cf = ceil(faces(i,j)/height);
                rf = mod(faces(i,j)-1,height)+1;
                if mask(rf,cf)
                    neighbor_face_flag(i) = false;
                    break;
                end
            end
        end

        if sum(neighbor_face_flag) == 0
            vertex_flag(vertex_index) = false;
            continue;
        end

        % guarantee that there exist at least a horizontal and a vertical
            % edges in the faces
        % alerady guaranteed by morphlogical operation of the mask
%         valid_faces = faces(neighbor_face_flag,:);
%         vflag = false;
%         hflag = false;
%         for i = 1 : size(valid_faces,1)
%             for j = 2 :3
%                 cf = ceil(valid_faces(i,j)/height);
%                 rf = mod(valid_faces(i,j)-1,height)+1;
%                 if cf == c
%                     vflag = true; % there exists a vertical edge
%                 end
%                 if rf == r
%                     hflag = true; % there exists a horizontal edge
%                 end
%                 if vflag && hflag
%                     break;
%                 end
%             end
%             if vflag && hflag
%                 break;
%             end
%         end
% 
%         if ~(vflag && hflag)
%             vertex_flag(vertex_index) = false;
%             face_flag(face_indices) = false;
%             continue;
%         end
        
        % find existing edges
        edge_indices_e = edge_indices(neighbor_face_flag,:);
        edge_flag(edge_indices_e(:)) = true;
        % store the face indices
        fov_index_map(r,c,1:sum(neighbor_face_flag)) = face_indices(neighbor_face_flag);
    end
end

% update index map
num_valid_vertex = sum(vertex_flag);
vertex_indices = zeros(height*width,1);
vertex_indices(vertex_flag) = 1 : num_valid_vertex;
vertex_index_map = reshape(vertex_indices, height, width);
V_INDEX_MAP = vertex_index_map;
%
vertices = reshape(P3MAP, [], 3); % (1,1) -> 1; (2,1) -> 2
vertices = vertices(vertex_flag,:);
VERTICES = vertices;
% update vertex indices in the faces
% make sure the first index belongs to the right angle
odd_faces = zeros(height-1,width-1,3);
odd_faces(:,:,1) = vertex_index_map(1:height-1,1:width-1);
odd_faces(:,:,2) = vertex_index_map(2:height,1:width-1);
odd_faces(:,:,3) = vertex_index_map(1:height-1,2:width);
odd_faces = reshape(odd_faces, [], 3);
even_faces = zeros(height-1,width-1,3);
even_faces(:,:,3) = vertex_index_map(2:height,1:width-1);
even_faces(:,:,1) = vertex_index_map(2:height,2:width);
even_faces(:,:,2) = vertex_index_map(1:height-1,2:width);
even_faces = reshape(even_faces, [], 3);
faces = zeros(size(odd_faces,1)+size(even_faces,1),3);
faces(1:2:end,:) = odd_faces;
faces(2:2:end,:) = even_faces;
%
faces = faces(face_flag,:);
FACES = faces;
% update the indices in fov_index_map
face_index_update = zeros(length(face_flag),1);
face_index_update(face_flag) = 1 : sum(face_flag);
fov_index_array = fov_index_map(:);
fov_index_array(fov_index_array~=0) = face_index_update(fov_index_array(fov_index_array~=0));
fov_index_map(:) = fov_index_array;
FOV_INDEX_MAP = fov_index_map;
% update vertex indices in the edges
edges_h = zeros(height,(width-1),2);
edges_h(:,:,1) = vertex_index_map(:,1:width-1);
edges_h(:,:,2) = vertex_index_map(:,2:width);
edges_h = reshape(edges_h, [], 2);
edges_v = zeros((height-1),width,2);
edges_v(:,:,1) = vertex_index_map(1:height-1,:);
edges_v(:,:,2) = vertex_index_map(2:height,:);
edges_v = reshape(edges_v, [], 2);
edges_d = zeros((height-1),(width-1),2);
edges_d(:,:,1) = vertex_index_map(2:height,1:width-1);
edges_d(:,:,2) = vertex_index_map(1:height-1,2:width);
edges_d = reshape(edges_d, [], 2);
% edges = zeros(size(edges_h,1)+size(edges_v,1)+size(edges_d,1),2);
% edges(1:size(edges_h,1),:) = edges_h;
% edges(size(edges_h,1)+1:size(edges_h,1)+size(edges_v,1),:) = edges_v;
% edges(size(edges_h,1)+size(edges_v,1)+1:end,:) = edges_d;
%
HEDGES = edges_h(edge_flag(1:size(edges_h,1)),:);
VEDGES = edges_v(edge_flag(size(edges_h,1)+1:size(edges_h,1)+size(edges_v,1)),:);
DEDGES = edges_d(edge_flag(size(edges_h,1)+size(edges_v,1)+1:end),:);
%
% edges = edges(edge_flag,:);
% EDGES = edges;

if CF
    colors = double(reshape(CMAP, [], 3))/255;
    colors = colors(vertex_flag,:);
    COLORS = colors;
else
    COLORS = [];
end

% MESH = surfaceMesh(vertices,faces,VertexColors=colors);

end
%%
function [edge_indices] = findNeighborEdges(r,c,h,w)

num_h = (w-1)*h; % number of horizontal edges
num_v = w*(h-1);
% num_d = (w-1)*(h-1);

if r < 1 || c < 1 || r > h || c > w
    edge_indices = [];
    disp("error: out of boundary");
elseif r == 1 && c == 1
    edge_indices = [(c-1)*h+r;...
        num_h+(c-1)*h+r];
elseif r == h && c == w
    edge_indices = [(c-1-1)*h+r;...
        num_h+(c-1)*(h-1)+(r-1)];
elseif r == h && c == 1
    edge_indices = [(c-1)*h+r;...
        num_h+(c-1)*(h-1)+(r-1);...
        num_h+num_v+(c-1)*(h-1)+(r-1)];
elseif r == 1 && c == w
    edge_indices = [(c-1-1)*h+r;
        num_h+(c-1)*(h-1)+r;
        num_h+num_v+(c-1-1)*(h-1)+(r+1-1)];
elseif r == 1
    edge_indices = [(c-1-1)*h+r;
        (c-1)*h+r;
        num_h+(c-1)*(h-1)+r;
        num_h+num_v+(c-1-1)*(h-1)+(r+1-1)];
elseif r == h
    edge_indices = [(c-1-1)*h+r;...
        num_h+(c-1)*(h-1)+(r-1);...
        num_h+num_v+(c-1)*(h-1)+(r-1)];
elseif c == 1
    edge_indices = [(c-1)*h+r;...
        num_h+(c-1)*(h-1)+r;...
        num_h+num_v+(c-1)*(h-1)+(r-1)];
elseif c == w
    edge_indices = [(c-1-1)*h+r;...
        num_h+(c-1)*(h-1)+(r-1);...
        num_h+num_v+(c-1-1)*(h-1)+(r+1-1)];
else
    edge_indices = [(c-1-1)*h+r;...
        (c-1)*h+r;...
        num_h+(c-1)*(h-1)+(r-1);...
        num_h+(c-1)*(h-1)+r;...
        num_h+num_v+(c-1-1)*(h-1)+(r+1-1);...
        num_h+num_v+(c-1)*(h-1)+(r-1)];
end

end
% %%
% function [index] = neighborEdgeIndexGenerator(r,c,i)
% 
% %      |    /
% %      3  5
% %      |/
% % --1--o--2--
% %     /|
% %   6  4
% % /    |
% if i == 1
%     index = (c-1-1)*h+r;
% else
%     index = [];
%     disp("error: wrong index");
% end
% 
% end
%%
function [face_indices,faces] = findNeighborFaces(r,c,h,w)

if r < 1 || c < 1 || r > h || c > w
    face_indices = [];
    faces = [];
    disp("error: out of boundary");
elseif r == 1 && c == 1
    face_indices = (c-1)*2*(h-1)+2*r-1;
    faces = [(c-1)*h+r,(c-1)*h+r+1,(c-1+1)*h+r];
elseif r == h && c == w
    face_indices = (c-1-1)*2*(h-1)+2*(r-1);
    faces = [(c-1)*h+r,(c-1)*h+r-1,(c-1-1)*h+r];
elseif r == h && c == 1
    face_indices = [(c-1)*2*(h-1)+2*(r-1)-1;...
        (c-1)*2*(h-1)+2*(r-1)];
    faces = [(c-1)*h+r,(c+1-1)*h+(r-1),(c-1)*h+(r-1);...
        (c-1)*h+r,(c+1-1)*h+r,(c+1-1)*h+(r-1)];
elseif r == 1 && c == w
    face_indices = [(c-1-1)*2*(h-1)+2*r-1;...
        (c-1-1)*2*(h-1)+2*r];
    faces = [(c-1)*h+r,(c-1-1)*h+r,(c-1-1)*h+r+1;...
        (c-1)*h+r,(c-1-1)*h+r+1,(c-1)*h+r+1];
elseif r == 1
    face_indices = [(c-1-1)*2*(h-1)+2*r-1;...
        (c-1-1)*2*(h-1)+2*r;...
        (c-1)*2*(h-1)+2*r-1];
    faces = [(c-1)*h+r,(c-1-1)*h+r,(c-1-1)*h+r+1;...
        (c-1)*h+r,(c-1-1)*h+r+1,(c-1)*h+r+1;...
        (c-1)*h+r,(c-1)*h+r+1,(c+1-1)*h+r];
elseif r == h
    face_indices = [(c-1-1)*2*(h-1)+2*(r-1);...
        (c-1)*2*(h-1)+2*(r-1)-1;...
        (c-1)*2*(h-1)+2*(r-1)];
    faces = [(c-1)*h+r,(c-1)*h+(r-1),(c-1-1)*h+r;...
        (c-1)*h+r,(c+1-1)*h+(r-1),(c-1)*h+(r-1);...
        (c-1)*h+r,(c+1-1)*h+r,(c+1-1)*h+(r-1)];
elseif c == 1
    face_indices = [(c-1)*2*(h-1)+2*(r-1)-1;...
        (c-1)*2*(h-1)+2*(r-1);...
        (c-1)*2*(h-1)+2*r-1];
    faces = [(c-1)*h+r,(c+1-1)*h+(r-1),(c-1)*h+(r-1);...
        (c-1)*h+r,(c+1-1)*h+r,(c+1-1)*h+(r-1);...
        (c-1)*h+r,(c-1)*h+r+1,(c+1-1)*h+r];
elseif c == w
    face_indices = [(c-1-1)*2*(h-1)+2*(r-1);...
        (c-1-1)*2*(h-1)+2*r-1;...
        (c-1-1)*2*(h-1)+2*r];
    faces = [(c-1)*h+r,(c-1)*h+(r-1),(c-1-1)*h+r;...
        (c-1)*h+r,(c-1-1)*h+r,(c-1-1)*h+r+1;...
        (c-1)*h+r,(c-1-1)*h+r+1,(c-1)*h+r+1];
else
    face_indices = [(c-1-1)*2*(h-1)+2*(r-1);...
        (c-1-1)*2*(h-1)+2*r-1;...
        (c-1-1)*2*(h-1)+2*r;...
        (c-1)*2*(h-1)+2*(r-1)-1;...
        (c-1)*2*(h-1)+2*(r-1);...
        (c-1)*2*(h-1)+2*r-1];
    faces = [(c-1)*h+r,(c-1)*h+(r-1),(c-1-1)*h+r;...
        (c-1)*h+r,(c-1-1)*h+r,(c-1-1)*h+r+1;...
        (c-1)*h+r,(c-1-1)*h+r+1,(c-1)*h+r+1;...
        (c-1)*h+r,(c+1-1)*h+(r-1),(c-1)*h+(r-1);...
        (c-1)*h+r,(c+1-1)*h+r,(c+1-1)*h+(r-1);...
        (c-1)*h+r,(c-1)*h+r+1,(c+1-1)*h+r];
end

end
%%
function [face_indices,faces,edge_indices] = findNeighborFacesEdges(r,c,h,w)

num_h = (w-1)*h; % number of horizontal edges
num_v = w*(h-1);

if r < 1 || c < 1 || r > h || c > w
    face_indices = [];
    faces = [];
    edge_indices = [];
    disp("error: out of boundary");
elseif r == 1 && c == 1
    face_indices = (c-1)*2*(h-1)+2*r-1;
    faces = [(c-1)*h+r,(c-1)*h+r+1,(c-1+1)*h+r];
    edge_indices = [(c-1)*h+r,num_h+(c-1)*(h-1)+r];
elseif r == h && c == w
    face_indices = (c-1-1)*2*(h-1)+2*(r-1);
    faces = [(c-1)*h+r,(c-1)*h+r-1,(c-1-1)*h+r];
    edge_indices = [(c-1-1)*h+r,num_h+(c-1)*(h-1)+(r-1)];
elseif r == h && c == 1
    face_indices = [(c-1)*2*(h-1)+2*(r-1)-1;...
        (c-1)*2*(h-1)+2*(r-1)];
    faces = [(c-1)*h+r,(c+1-1)*h+(r-1),(c-1)*h+(r-1);...
        (c-1)*h+r,(c+1-1)*h+r,(c+1-1)*h+(r-1)];
    edge_indices = [num_h+(c-1)*(h-1)+(r-1),num_h+num_v+(c-1)*(h-1)+(r-1);...
        num_h+num_v+(c-1)*(h-1)+(r-1),(c-1)*h+r];
elseif r == 1 && c == w
    face_indices = [(c-1-1)*2*(h-1)+2*r-1;...
        (c-1-1)*2*(h-1)+2*r];
    faces = [(c-1)*h+r,(c-1-1)*h+r,(c-1-1)*h+r+1;...
        (c-1)*h+r,(c-1-1)*h+r+1,(c-1)*h+r+1];
    edge_indices = [(c-1-1)*h+r,num_h+num_v+(c-1-1)*(h-1)+(r+1-1);...
        num_h+num_v+(c-1-1)*(h-1)+(r+1-1),num_h+(c-1)*(h-1)+r];
elseif r == 1
    face_indices = [(c-1-1)*2*(h-1)+2*r-1;...
        (c-1-1)*2*(h-1)+2*r;...
        (c-1)*2*(h-1)+2*r-1];
    faces = [(c-1)*h+r,(c-1-1)*h+r,(c-1-1)*h+r+1;...
        (c-1)*h+r,(c-1-1)*h+r+1,(c-1)*h+r+1;...
        (c-1)*h+r,(c-1)*h+r+1,(c+1-1)*h+r];
    edge_indices = [(c-1-1)*h+r,num_h+num_v+(c-1-1)*(h-1)+(r+1-1);...
        num_h+num_v+(c-1-1)*(h-1)+(r+1-1),num_h+(c-1)*(h-1)+r;...
       (c-1)*h+r,num_h+(c-1)*(h-1)+r];
elseif r == h
    face_indices = [(c-1-1)*2*(h-1)+2*(r-1);...
        (c-1)*2*(h-1)+2*(r-1)-1;...
        (c-1)*2*(h-1)+2*(r-1)];
    faces = [(c-1)*h+r,(c-1)*h+(r-1),(c-1-1)*h+r;...
        (c-1)*h+r,(c+1-1)*h+(r-1),(c-1)*h+(r-1);...
        (c-1)*h+r,(c+1-1)*h+r,(c+1-1)*h+(r-1)];
    edge_indices = [(c-1-1)*h+r,num_h+(c-1)*(h-1)+(r-1);...
        num_h+(c-1)*(h-1)+(r-1),num_h+num_v+(c-1)*(h-1)+(r-1);...
        num_h+num_v+(c-1)*(h-1)+(r-1),(c-1)*h+r];
elseif c == 1
    face_indices = [(c-1)*2*(h-1)+2*(r-1)-1;...
        (c-1)*2*(h-1)+2*(r-1);...
        (c-1)*2*(h-1)+2*r-1];
    faces = [(c-1)*h+r,(c+1-1)*h+(r-1),(c-1)*h+(r-1);...
        (c-1)*h+r,(c+1-1)*h+r,(c+1-1)*h+(r-1);...
        (c-1)*h+r,(c-1)*h+r+1,(c+1-1)*h+r];
    edge_indices = [num_h+(c-1)*(h-1)+(r-1),num_h+num_v+(c-1)*(h-1)+(r-1);...
        num_h+num_v+(c-1)*(h-1)+(r-1),(c-1)*h+r;...
        (c-1)*h+r,num_h+(c-1)*(h-1)+r];
elseif c == w
    face_indices = [(c-1-1)*2*(h-1)+2*(r-1);...
        (c-1-1)*2*(h-1)+2*r-1;...
        (c-1-1)*2*(h-1)+2*r];
    faces = [(c-1)*h+r,(c-1)*h+(r-1),(c-1-1)*h+r;...
        (c-1)*h+r,(c-1-1)*h+r,(c-1-1)*h+r+1;...
        (c-1)*h+r,(c-1-1)*h+r+1,(c-1)*h+r+1];
    edge_indices = [(c-1-1)*h+r,num_h+(c-1)*(h-1)+(r-1);...
        (c-1-1)*h+r,num_h+num_v+(c-1-1)*(h-1)+(r+1-1);...
        num_h+num_v+(c-1-1)*(h-1)+(r+1-1),num_h+(c-1)*(h-1)+r];
else
    face_indices = [(c-1-1)*2*(h-1)+2*(r-1);...
        (c-1-1)*2*(h-1)+2*r-1;...
        (c-1-1)*2*(h-1)+2*r;...
        (c-1)*2*(h-1)+2*(r-1)-1;...
        (c-1)*2*(h-1)+2*(r-1);...
        (c-1)*2*(h-1)+2*r-1];
    faces = [(c-1)*h+r,(c-1)*h+(r-1),(c-1-1)*h+r;...
        (c-1)*h+r,(c-1-1)*h+r,(c-1-1)*h+r+1;...
        (c-1)*h+r,(c-1-1)*h+r+1,(c-1)*h+r+1;...
        (c-1)*h+r,(c+1-1)*h+(r-1),(c-1)*h+(r-1);...
        (c-1)*h+r,(c+1-1)*h+r,(c+1-1)*h+(r-1);...
        (c-1)*h+r,(c-1)*h+r+1,(c+1-1)*h+r];
    edge_indices = [(c-1-1)*h+r,num_h+(c-1)*(h-1)+(r-1);...
        (c-1-1)*h+r,num_h+num_v+(c-1-1)*(h-1)+(r+1-1);...
        num_h+num_v+(c-1-1)*(h-1)+(r+1-1),num_h+(c-1)*(h-1)+r;...
        num_h+(c-1)*(h-1)+(r-1),num_h+num_v+(c-1)*(h-1)+(r-1);...
        num_h+num_v+(c-1)*(h-1)+(r-1),(c-1)*h+r;...
        (c-1)*h+r,num_h+(c-1)*(h-1)+r];
end

end
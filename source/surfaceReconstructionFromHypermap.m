% CF : color flag

%%
function [MESH] = surfaceReconstructionFromHypermap(HYPERMAP,CF)

P3MAP = HYPERMAP.point_3d_map;
MASK = logical(~HYPERMAP.vertex_index_map);
CMAP = HYPERMAP.color_map;

height = size(P3MAP,1);
width = size(P3MAP,2);

vertex_flag = true(width*height,1);
face_flag = true(2*(width-1)*(height-1),1);
for c = 1 : width
    for r = 1 : height
        vertex_index = (c-1)*height+r;
        [face_indices,~,~] = findNeighborFacesEdges(r,c,height,width);
        % check mask
        if MASK(r,c)
            vertex_flag(vertex_index) = false;
            face_flag(face_indices) = false;
            continue;
        end
    end
end

% update index map
num_valid_vertex = sum(vertex_flag);
vertex_indices = zeros(height*width,1);
vertex_indices(vertex_flag) = 1 : num_valid_vertex;
vertex_index_map = reshape(vertex_indices, height, width);
%
vertices = reshape(P3MAP, [], 3); % (1,1) -> 1; (2,1) -> 2
vertices = vertices(vertex_flag,:);
% update vertex indices in the faces
odd_faces = zeros(height-1,width-1,3);
odd_faces(:,:,1) = vertex_index_map(1:height-1,1:width-1);
odd_faces(:,:,2) = vertex_index_map(2:height,1:width-1);
odd_faces(:,:,3) = vertex_index_map(1:height-1,2:width);
odd_faces = reshape(odd_faces, [], 3);
even_faces = zeros(height-1,width-1,3);
even_faces(:,:,1) = vertex_index_map(2:height,1:width-1);
even_faces(:,:,2) = vertex_index_map(2:height,2:width);
even_faces(:,:,3) = vertex_index_map(1:height-1,2:width);
even_faces = reshape(even_faces, [], 3);
faces = zeros(size(odd_faces,1)+size(even_faces,1),3);
faces(1:2:end,:) = odd_faces;
faces(2:2:end,:) = even_faces;
%
faces = faces(face_flag,:);

if CF
    colors = double(reshape(CMAP, [], 3))/255;
    colors = colors(vertex_flag,:);
    MESH = surfaceMesh(vertices,faces,VertexColors=colors);
else
    MESH = triangulation(faces,vertices);
end

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
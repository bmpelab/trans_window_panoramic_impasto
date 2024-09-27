% BSFS : backward scene flow vectors
% BINDICES : == n means n-th boundary
% PL : left projection matrix
% SR : sampling rate
% OFFSET : image offset

% P3MAP : point 3d map
% CMAP : color map, for color, closest interpolation
% DMAP : derivative map
% BSFMAP : backward scene flow map
% MASK : == true means undefined
% BIMAP : == n means n-th boundary

%%
function [P3MAP,DMAP,CMAP,BSFMAP,BOFMAP,MASK,BIMAP] = parameterization_ver7(H,W,HC,WC,VERTICES,FACES,DERIVATIVES,COLORS,BSFS,BOFS,BINDICES,PL,SR,OFFSET)

dH = round((HC-H)/2);
dW = round((WC-W)/2);

F = PL(1,1);
CX = PL(1,3);
CY = PL(2,3);

VN = vertexNormal(triangulation(FACES,VERTICES));
VV = -VERTICES./sqrt(sum(VERTICES.^2,2));
front_flags = sum(VV.*VN,2)>0;

% project each face onto the canonical space
    % each face forms a bounding box
    % keep the face id within the bounding box
distance_map = inf(HC,WC);
cmap = uint8(zeros(HC,WC,3));
p3map = zeros(HC,WC,3);
dmap = zeros(HC,WC,6);
bsfmap = zeros(HC,WC,3);
bofmap = zeros(HC,WC,2);
mask = true(HC,WC);
bimap = zeros(HC,WC);
for i = 1 : size(FACES,1)
    % the first index belongs to the right angle
    face = FACES(i,:);
    front_flag = front_flags(face);
    if sum(front_flag) ~= 3
        continue;
    end
    vertices = VERTICES(face,:);
    % if the angle between the face normal and the face center vector is larger than 90
    pixels = PL*[vertices';ones(1,3)];
    pixels = pixels./pixels(3,:);
    pixels(1,:) = pixels(1,:)-OFFSET;
    pixels = pixels/SR;
    pixels = pixels(1:2,:)';
    % convert pixels into the canonical coordinates
    % (y=0,x=0) -> (r=1,c=1) -> round(height/2)+1,round(width/2)+1
    % (y=height-1,x=width-1) -> (r=height,c=width) -> round(height/2)+height,round(width/2)+width
    pixels_canonical = zeros(size(pixels));
    pixels_canonical(:,1) = pixels(:,1)+dW; % x
    pixels_canonical(:,2) = pixels(:,2)+dH; % y
    % boundary box
    ra = ceil(min(pixels_canonical(:,2))+1);
    rb = floor(max(pixels_canonical(:,2))+1);
    ca = ceil(min(pixels_canonical(:,1))+1);
    cb = floor(max(pixels_canonical(:,1))+1);
    %
    if ra < 1
        ra = 1;
    elseif ra > HC
        continue;
    end
    if rb < 1
        continue;
    elseif rb > HC
        rb = HC;
    end
    if ca < 1
        ca = 1;
    elseif ca > WC
        continue;
    end
    if cb < 1
        continue;
    elseif cb > WC
        cb = WC;
    end
    %
    if ra > rb
        continue;
    end
    if ca > cb
        continue;
    end
    %
    for r = ra : rb
        for c = ca : cb
            % within the triangle formed by the projection of the face
            if isInsideTriangle([c-1,r-1,0], [pixels_canonical(1,:),0], [pixels_canonical(2,:),0], [pixels_canonical(3,:),0])
                % calculate the intersection
                % coordinates in image coordinates
                u = ((c-1)-dW)*SR+OFFSET;
                v = ((r-1)-dH)*SR;
                line_points = [(u-CX)/F,(v-CY)/F,1;0,0,0];
                [intersection_point] = lineFaceIntersection(vertices,line_points,false);
                distance_to_camera = norm(intersection_point);
                if distance_to_camera < distance_map(r,c) % do interpolation
                    distance_map(r,c) = distance_to_camera;
                    % calculate the weights for interpolation using the distances
                    distances = sqrt(sum((vertices-intersection_point).^2,2));
                    % nearest-neighbor interpolation for color
                    [min_d,min_i] = min(distances);
                    cmap(r,c,:) = COLORS(face(min_i),:);
                    % do not interpolate the point 3d
                    % just use the intersection point
                    p3map(r,c,:) = intersection_point;
                    if min_d ~= 0
                        % inverse-distance interpolation
                        inverse_distances = 1./distances;
                        weights = inverse_distances/sum(inverse_distances);
                        weights = weights(:)'; % 1 by 3
                        %
                        derivatives = DERIVATIVES(face,:); % 3 by 6
                        derivatives_ = zeros(size(derivatives));
                        for j = 1 : 3
                            % calculate the intersection
                            % coordinates in image coordinates
                            line_point = [0,0,0];
                            % to the right (horizontal)
                            u = ((c-1+1)-dW)*SR+OFFSET;
                            v = ((r-1)-dH)*SR;
                            line_point2 = [(u-CX)/F,(v-CY)/F,1];
                            line_direction = line_point2/norm(line_point2);
                            plane_point = intersection_point;
                            plane_normal = cross(derivatives(j,1:3),derivatives(j,4:6));
                            plane_normal = plane_normal/norm(plane_normal);
                            [intersection_point2] = linePlaneIntersection(line_direction,line_point,plane_normal,plane_point);
                            derivatives_(j,1:3) = intersection_point2 - intersection_point;
                            % to the bottom (vertical)
                            u = ((c-1)-dW)*SR+OFFSET;
                            v = ((r-1+1)-dH)*SR;
                            line_point2 = [(u-CX)/F,(v-CY)/F,1];
                            line_direction = line_point2/norm(line_point2);
                            [intersection_point2] = linePlaneIntersection(line_direction,line_point,plane_normal,plane_point);
                            derivatives_(j,4:6) = intersection_point2 - intersection_point;
                        end
                        dmap(r,c,:) = [weights*derivatives_(:,1:3),weights*derivatives_(:,4:6)];
                        %
                        bsfs = BSFS(face,:);
                        bsfmap(r,c,:) = weights*bsfs;
                        %
                        bofs = BOFS(face,:);
                        bofmap(r,c,:) = weights*bofs;
                    else
                        %
                        derivatives = DERIVATIVES(face(min_i),:); % 1 by 6
                        derivatives_ = zeros(size(derivatives));
                        for j = 1 : 1
                            % calculate the intersection
                            % coordinates in image coordinates
                            line_point = [0,0,0];
                            % to the right (horizontal)
                            u = ((c-1+1)-dW)*SR+OFFSET;
                            v = ((r-1)-dH)*SR;
                            line_point2 = [(u-CX)/F,(v-CY)/F,1];
                            line_direction = line_point2/norm(line_point2);
                            plane_point = intersection_point;
                            plane_normal = cross(derivatives(j,1:3),derivatives(j,4:6));
                            plane_normal = plane_normal/norm(plane_normal);
                            [intersection_point2] = linePlaneIntersection(line_direction,line_point,plane_normal,plane_point);
                            derivatives_(j,1:3) = intersection_point2 - intersection_point;
                            % to the bottom (vertical)
                            u = ((c-1)-dW)*SR+OFFSET;
                            v = ((r-1+1)-dH)*SR;
                            line_point2 = [(u-CX)/F,(v-CY)/F,1];
                            line_direction = line_point2/norm(line_point2);
                            [intersection_point2] = linePlaneIntersection(line_direction,line_point,plane_normal,plane_point);
                            derivatives_(j,4:6) = intersection_point2 - intersection_point;
                        end
                        dmap(r,c,:) = derivatives_;
                        bsfmap(r,c,:) = BSFS(face(min_i),:);
                        bofmap(r,c,:) = BOFS(face(min_i),:);
                    end
                    mask(r,c) = false;
                end
            end
        end
    end
    %
    for j = 1 : 3
        if BINDICES(face(j))>inf
            continue;
        end
        %
        r = round(pixels_canonical(j,2)+1);
        c = round(pixels_canonical(j,1)+1);
        %
        if r < 1 || r > HC || c < 1 || c > WC
            continue;
        end
        % skip if there exists info
        if ~mask(r,c)
            continue;
        end
        %
        vertex = vertices(j,:);
        p3map(r,c,:) = vertex;
        cmap(r,c,:) = COLORS(face(j),:); % [0 255 0];
        %
        derivatives = DERIVATIVES(face(j),:); % 1 by 6
        derivatives_ = zeros(size(derivatives));
        % calculate the intersection
        % coordinates in image coordinates
        line_point = [0,0,0];
        % to the right (horizontal)
        u = ((c-1+1)-dW)*SR+OFFSET;
        v = ((r-1)-dH)*SR;
        line_point2 = [(u-CX)/F,(v-CY)/F,1];
        line_direction = line_point2/norm(line_point2);
        plane_point = vertex;
        plane_normal = cross(derivatives(1:3),derivatives(4:6));
        plane_normal = plane_normal/norm(plane_normal);
        [intersection_point2] = linePlaneIntersection(line_direction,line_point,plane_normal,plane_point);
        derivatives_(1:3) = intersection_point2 - vertex;
        % to the bottom (vertical)
        u = ((c-1)-dW)*SR+OFFSET;
        v = ((r-1+1)-dH)*SR;
        line_point2 = [(u-CX)/F,(v-CY)/F,1];
        line_direction = line_point2/norm(line_point2);
        [intersection_point2] = linePlaneIntersection(line_direction,line_point,plane_normal,plane_point);
        derivatives_(4:6) = intersection_point2 - vertex;
        %
        dmap(r,c,:) = derivatives_;
        bsfmap(r,c,:) = BSFS(face(j),:);
        bofmap(r,c,:) = BOFS(face(j),:);
        mask(r,c) = false;
        bimap(r,c) = BINDICES(face(j))+1;
    end
end

% 
% % filter based on color map
% % clear the region where color diffused
% cmap_archives = cmap;
% for r = 3 : HC-2
%     for c = 3 : WC-2
%         if mask(r,c)
%             continue;
%         end
%         cmap_ = cmap_archives(r-2:r+2,c-2:c+2,:);
%         color = squeeze(cmap_(3,3,:))';
%         colors = reshape(cmap_,[],3);
%         if sum(sum(colors==color,2)==3) > 16
%             p3map(r,c,:) = 0;
%             cmap(r,c,:) = uint8(0);
%             dmap(r,c,:) = 0;
%             bsfmap(r,c,:) = 0;
%             mask(r,c) = true;
%         end
%     end
% end

BSFMAP = bsfmap;
BOFMAP = bofmap;
CMAP = cmap;
P3MAP = p3map;
DMAP = dmap;
MASK = mask;
BIMAP = bimap;

end
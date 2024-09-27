function [intersection_point] = lineFaceIntersection(face_points,line_points,checkflag)

% Extracting face points
A = face_points(1,:);
B = face_points(2,:);
C = face_points(3,:);

% Extracting line points
P = line_points(1,:);
Q = line_points(2,:);

% Vectors representing edges of the face
AB = B - A;
AC = C - A;

% Normal vector of the face
normal = cross(AB, AC);

% Equation of the plane of the face: Ax + By + Cz + D = 0
D = -dot(normal, A);

% Direction vector of the line
direction = Q - P;

% Parameter t for the line equation P + t*direction
t = -(dot(normal, P) + D) / dot(normal, direction);

% Intersection point
intersection_point = P + t*direction;

% Check if the intersection point lies within the triangle ABC
if checkflag
    if ~isInsideTriangle(intersection_point, A, B, C)
        intersection_point = []; % No intersection
    end
end

% figure;hold on;
% plot3([A(1),B(1),C(1),A(1)],[A(2),B(2),C(2),A(2)],[A(3),B(3),C(3),A(3)],'r-');
% plot3([P(1),Q(1)],[P(2),Q(2)],[P(3),Q(3)],'b-');
% plot3([intersection_point(1)],[intersection_point(2)],[intersection_point(3)],'k*');
end
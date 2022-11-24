function [intersection] = TriangleSegmentIntersection(seg, tri, epsilon)
% Compute whether there is an intersection between triangle and line segment.
% Input:
%   segment: 2*3 matrix. Each row is an end point of the line segment.
%   tri: 3*3 matrix. Each row is an vertex of the triangle.
%   epsilon: error tolerance.
% Ouput:
%   intersection: whether there is an intersection between triangle and 
%   line segment.
intersection = 0;
e1 = tri(1,:) - tri(2,:);
e2 = tri(1,:) - tri(3,:);
n = cross(e1,e2);
n = n/norm(n);

[int,check]=plane_line_intersect(n,tri(1,:),seg(1,:),seg(2,:));

if check == 1 % unique intersection in seg
    if PointInTriangle(int,tri,epsilon) > 0
        intersection = 1;
    end
end

if check == 2 % seg is in the same plane as tri
    if PointInTriangle(seg(1,:),tri,epsilon) || PointInTriangle(seg(2,:),tri,epsilon)
        intersection = 1;
    end
end

if check == 3 % unique intersection outside seg
    t1 = (seg(1,:)-int)*n'; % signed distance between seg(1,:) and plane containing tri.
    p1 = seg(1,:)-t1*n; % projection of seg(1,:) to triangle plane
    t1 = abs(t1);
    t2 = (seg(2,:)-int)*n';
    p2 = seg(2,:)-t2*n;
    t2 = abs(t2);
    if (t1<epsilon && PointInTriangle(p1,tri,epsilon) == 1)...
            || (t2 <epsilon && PointInTriangle(p2,tri,epsilon) == 1)
        intersection = 1;
    end
end
end
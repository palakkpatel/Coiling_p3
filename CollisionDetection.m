function [ret, tri_idx] = CollisionDetection(seg,V,F,epsilon)
% Check whether line segment collides with all the triangles.
% Input:
%   seg: 2*3 matrix. Each row is an end point of the line segment.
%   V: n*3 matrix. Each row is a vertex.
%   F: m*3 matrix. Each row is the indices of the vertices in V.
% Output:
%   ret: whether line segment collides with all the triangles.
%   tri: the triangle that collides with the line segment.
ret = 0;
n_tri = size(F,1);
tri_idx = 0;
dis = norm(seg(1,:)-seg(2,:))+epsilon;

for ii = 1:n_tri
    tri = [V(F(ii,1),:);V(F(ii,2),:);V(F(ii,3),:)];
    if norm(tri(1,:)-seg(1,:)) > 2*dis && norm(tri(2,:)-seg(1,:))>2*dis ...
            && norm(tri(3,:)-seg(1,:))>2*dis
        continue;
    end
    ret = TriangleSegmentIntersection(seg, tri, epsilon);
    if ret > 0
        tri_idx = ii;
        break;
    end
end

end
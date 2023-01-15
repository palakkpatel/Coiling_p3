function [new_seg] = SegmentRotation(seg,p,V,F,tri_idx,radius,epsilon)
% Rotate the line segment such that the segment inside the dome.
% Let's assume the second point is inside.
% Input:
%   seg: 2*3 matrix. Each row is an end point.
%   p: 1*3 vector. A reference point.
%   V: n*3 matrix. Each row is a vertex.
%   F: m*3 matrix. Each row is a triangle.
%   tri_idx: integer. The index of the triangle.
%   epsilon: distance needed between the segment and the triangle.
% Output:
%   new_seg: 2*3 matrix. Each row is an end point.

tri = [V(F(tri_idx,1),:);V(F(tri_idx,2),:);V(F(tri_idx,3),:)];
n1 = cross(tri(1,:)-tri(2,:),tri(1,:)-tri(3,:));
n2 = cross(seg(1,:)-seg(2,:),seg(1,:)-p);
n1 = n1/norm(n1);
n2 = n2/norm(n2);

if (seg(2,:)-tri(1,:))*n1'<0
    n1 = -n1;
end

for ii = 1:3
    tri(ii,:) = tri(ii,:) + n1*(radius+epsilon);
end

con1 = n1*tri(1,:)';
con2 = n2*seg(2,:)';

ax = (n1(2)*n2(3) - n2(2)*n1(3))/(n1(1)*n2(2) - n1(2)*n2(1));
bx = (n2(2)*con1 - n1(2)*con2)/(n1(1)*n2(2) - n1(2)*n2(1));
ay = (n1(1)*n2(3) - n2(1)*n1(3))/(n2(1)*n1(2) - n1(1)*n2(2));
by = (n2(1)*con1 - n1(1)*con2)/(n2(1)*n1(2) - n1(1)*n2(2));

bx = bx - seg(2,1);
by = by - seg(2,2);
R = norm(seg(2,:) - seg(1,:));

a = ax^2 + ay^2 + 1;
b = 2*ax*bx + 2*ay*by - 2*seg(2,3);
c = bx^2 + by^2 + seg(2,3)^2 - R^2;

new_seg(1,3) = (-b + sqrt(b^2 - 4*a*c))/(2*a);
new_seg(2,3) = (-b - sqrt(b^2 - 4*a*c))/(2*a);
new_seg(:,1) = ax*new_seg(:,3) + bx + seg(2,1);
new_seg(:,2) = ay*new_seg(:,3) + by + seg(2,2);

if norm(new_seg(1,:) - seg(1,:)) > norm(new_seg(2,:) - seg(1,:))
    new_seg(1,:) = new_seg(2,:);
end

new_seg(2,:) = seg(2,:);

end
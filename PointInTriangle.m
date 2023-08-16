function [ret] = PointInTriangle(p,tri, radius)
%% check whether point p lies in triangle tri of mesh
% Input:
%   p: 1*3 vector.
%   tri: 3*3 matrix. Each row is a point.
% Output:
%   ret: 0-outside, 1-inside.
ret = 1;
%% v1v2
% k1 = ((tri(3,:)-tri(1,:))*(tri(2,:)-tri(1,:))')/((tri(2,:)-tri(1,:))*(tri(2,:)-tri(1,:))');
% k2 = ((p-tri(1,:))*(tri(2,:)-tri(1,:))')/((tri(2,:)-tri(1,:))*(tri(2,:)-tri(1,:))');
% if (tri(1,:)+k1*(tri(2,:)-tri(1,:)) - tri(3,:))*(tri(1,:)+k2*(tri(2,:)-tri(1,:)) - p)' < 0
%     ret = 0;
% end


% e1 = tri(3,:)-tri(1,:);
% e1 = e1/norm(e1);
% e2 = (p-tri(1,:))*e1'*e1;
% if e1*e2'<0 && norm(e2)>radius
%     ret = 0;
% end

n = tri(2,:)-tri(1,:);
n = n/norm(n);
e1 = tri(3,:)-tri(1,:);
e2 = p-tri(1,:);
p1 = e1*n';
p2 = e2*n';
p1 = tri(1,:)+p1*n;
p2 = tri(1,:)+p2*n;
e1 = tri(3,:)-p1;
e2 = p-p2;
if e1*e2'< 0 && norm(e2) > radius
    ret = 0;
end

f1 = e1*e2';

%% v1v3
% k1 = ((tri(2,:)-tri(1,:))*(tri(3,:)-tri(1,:))')/((tri(3,:)-tri(1,:))*(tri(3,:)-tri(1,:))');
% k2 = ((p-tri(1,:))*(tri(3,:)-tri(1,:))')/((tri(3,:)-tri(1,:))*(tri(3,:)-tri(1,:))');
% if (tri(1,:)+k1*(tri(3,:)-tri(1,:)) - tri(2,:))*(tri(1,:)+k2*(tri(3,:)-tri(1,:)) - p)' < 0
%     ret = 0;
% end

% e1 = tri(2,:)-tri(1,:);
% e1 = e1/norm(e1);
% e2 = (p-tri(1,:))*e1'*e1;
% if e1*e2'<0 && norm(e2)>radius
%     ret = 0;
% end

n = tri(3,:)-tri(1,:);
n = n/norm(n);
e1 = tri(2,:)-tri(1,:);
e2 = p-tri(1,:);
p1 = e1*n';
p2 = e2*n';
p1 = tri(1,:)+p1*n;
p2 = tri(1,:)+p2*n;
e1 = tri(2,:)-p1;
e2 = p-p2;
if e1*e2'<0 && norm(e2)>radius
    ret = 0;
end

f2 = e1*e2';
%% v2v3
% k1 = ((tri(1,:)-tri(2,:))*(tri(3,:)-tri(2,:))')/((tri(3,:)-tri(2,:))*(tri(3,:)-tri(2,:))');
% k2 = ((p-tri(2,:))*(tri(3,:)-tri(2,:))')/((tri(3,:)-tri(2,:))*(tri(3,:)-tri(2,:))');
% if (tri(2,:)+k1*(tri(3,:)-tri(2,:)) - tri(1,:))*(tri(2,:)+k2*(tri(3,:)-tri(2,:)) - p)' < 0
%     ret = 0;
% end

% e1 = tri(1,:)-tri(2,:);
% e1 = e1/norm(e1);
% e2 = (p-tri(2,:))*e1'*e1;
% if e1*e2'<0 && norm(e2)>radius
%     ret = 0;

n = tri(3,:)-tri(2,:);
n = n/norm(n);
e1 = tri(1,:)-tri(2,:);
e2 = p-tri(2,:);
p1 = e1*n';
p2 = e2*n';
p1 = tri(2,:)+p1*n;
p2 = tri(2,:)+p2*n;
e1 = tri(1,:)-p1;
e2 = p-p2;
if e1*e2'<0 && norm(e2)>radius
    ret = 0;
end

f3 = e1*e2';

if f1<0 || f2<0 || f3<0
    if norm(p-tri(1,:))>radius || norm(p-tri(2,:))>radius || norm(p-tri(3,:))>radius
        ret = 0;
    end
end


end

function [collision,seg_idx] = collisionSeg(seg,meshpoint,radius,grid_s,grid_m)
% Check whether segment 'seg' intersects with 'meshpoint' on the Sac
% Input:
%   seg: 2x3 matrix.
%   meshpoint: mx3 matrix.
%   radius: radius of the coil.
% Output:
%   collision: 1 - collision; 0 - no collision.
%   seg_idx: if collision, the intersecting segment in meshpoint is
%               meshpoint(seg_idx-1:seg_idx,:).

collision = 0;
for seg_idx = 2:(size(meshpoint,1)-1)
    if sum(grid_s(1,:) == grid_m(seg_idx,:)) < 3 && ...
        sum(grid_s(2,:) == grid_m(seg_idx,:)) < 3 && ...
        sum(grid_s(1,:) == grid_m(seg_idx-1,:))<3 && ...
        sum(grid_s(2,:) == grid_m(seg_idx-1,:))<3
%     if  sum(grid_s(2,:) == grid_m(seg_idx,:)) < 3 && ...
%         sum(grid_s(2,:) == grid_m(seg_idx-1,:))<3
        continue;
    end
    if norm(seg(2,:)-meshpoint(seg_idx,:))>4*radius
        continue;
    end
    if DistBetween2Segment(seg(1,:),seg(2,:),...
            meshpoint(seg_idx-1,:),meshpoint(seg_idx,:))<2.5*radius
        collision = 1;
        break;
    end
end

end
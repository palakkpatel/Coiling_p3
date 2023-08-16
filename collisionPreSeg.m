function [collision, seg_idx] = collisionPreSeg(seg,temp_p,temp_h,temp_c,radius,v_min,v_max,grid_n)
% Check collision with previous Coils
% Input:
%   seg: 2*3 matrix. Current segment.
%   temp_p: n*3 matrix. deployment of previous coil.
%   temp_h: (3n)*2 matrix. helping table of previous coil.
%   temp_c: v_cnt*3 matrix. count info of helping table.
%   radius: scalar. radius of segment.
% Ouput:    
% collision: scalar. 0: no collision. 1: collision.
% seg_idx: scalar. The index of segment collides with seg.

collision = 0;
n_seg = size(temp_p,1);
seg_idx = 0;
dis = norm(seg(1,:)-seg(2,:));
grid_s = [update_grid(seg(1,:),v_min,v_max,grid_n); ...
    update_grid(seg(2,:),v_min,v_max,grid_n)];


for ii = 1:n_seg-1
    grid_m = temp_p(ii:ii+1,4:6);
    if sum(grid_s(1,:) == grid_m(1,:)) < 3 && ...
        sum(grid_s(2,:) == grid_m(1,:)) < 3 && ...
        sum(grid_s(1,:) == grid_m(2,:))<3 && ...
        sum(grid_s(2,:) == grid_m(2,:))<3
        continue;
    end
    if norm(seg(1,:)-temp_p(ii,1:3)) > 4*radius
        continue;
    end
    if DistBetween2Segment(seg(1,:),seg(2,:),...
            temp_p(ii,1:3),temp_p(ii+1,1:3))<2.5*radius
        collision = 1;
        seg_idx = ii;
        break;
    end
end

end
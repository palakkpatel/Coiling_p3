function [help_t, cnt_t] = previousCoils(meshpoint,grid_location)
% Compute the helping table for previously deployed coils.
% This table helps when computing the collision between current coil and
% previous coils.
% Input:
%   meshpoint: n*3 matrix. The deployment of coil.
%   grid_location: n*3 matrix. The grid_location of each meshpoint.
% Output:
%   help_t: (3n)*2 matrix. The helping table for this coil.
%   cnt_t: v_cnt*1 vector, store the ending index of each grid index.

cnt_p = size(meshpoint,1);
maxgrid = max(grid_location);
grid_location(:,2) = grid_location(:,2) + maxgrid(1);
grid_location(:,3) = grid_location(:,3) + maxgrid(1) + maxgrid(2);

help_t = zeros(3*cnt_p,2);
help_t(1:cnt_p,1) = grid_location(:,1);
help_t(1:cnt_p,2) = [1:cnt_p]';
help_t(cnt_p+1:2*cnt_p,1) = grid_location(:,2);
help_t(cnt_p+1:2*cnt_p,2) = [1:cnt_p]';
help_t(2*cnt_p+1:3*cnt_p,1) = grid_location(:,3);
help_t(2*cnt_p+1:3*cnt_p,2) = [1:cnt_p]';

[temp, idx] = sort(help_t(:,1));
help_t(:,1) = temp;
help_t(:,2) = help_t(idx,2);
v_cnt = max(help_t(:,1));
cnt_t = zeros(v_cnt,1);
for ii = 1:v_cnt
    cnt_t(ii) = size(find(help_t(:,1) == ii),1);
end

for ii = 2:v_cnt
    cnt_t(ii) = cnt_t(ii-1)+cnt_t(ii);
end

end
function [tri_id] = tri_neighbour(tri_idx,fout,cnt, help_t, cnt_t)
% find the cnt neighbours of tri_idx.
% the size of tri_id should be no less than cnt.
% Input:
%   tri_idx: scalar, the index of target triangle.
%   fout: f_cnt*3 matrix, the triangle mesh of the dome.
%   cnt: scalar, the number of neighbours wanted.
%   help_t: helping table.
%   cnt_t: count of the vertices in helping table.
% Ouput:
%   tri_id: indeces of neighbours in fout.

p = fout(tri_idx,:);
idx = zeros(3,2);
for ii = 1:3
    if p(ii) == 1
        idx(ii,1) = 1;
    else 
        idx(ii,1) = cnt_t(p(ii)-1)+1;
    end
    idx(ii,2) = cnt_t(p(ii));
end

tri_id = [help_t(idx(1,1):idx(1,2),2);...
    help_t(idx(2,1):idx(2,2),2);...
    help_t(idx(3,1):idx(3,2),2)];
tri_id = unique(tri_id);
tri_cnt = size(tri_id,1);
temp = [];

while tri_cnt < cnt
    for ii = 1:tri_cnt
        p = fout(tri_id(ii),:);
        idx = zeros(3,2);
        for jj = 1:3
            if p(jj) == 1
                idx(jj,1) = 1;
            else 
                idx(jj,1) = cnt_t(p(jj)-1)+1;
            end
            idx(jj,2) = cnt_t(p(jj));
        end
        temp = [temp; help_t(idx(1,1):idx(1,2),2);...
                help_t(idx(2,1):idx(2,2),2);...
                help_t(idx(3,1):idx(3,2),2)];
        temp = unique(temp);
    end
    tri_id = [tri_id;temp];
    tri_id = unique(tri_id);
    tri_cnt = size(tri_id,1);
end
end
function [help_t,cnt_t] = neighbour_table(fout)
% build the helping table to find noughbours
% Input:
%   fout: f_cnt*3 matrix, triangle mesh of the dome.
% Output:
%   help_t: (3*f_cnt)*2 matrix, stores the vertex indeces of the triangle in
%   order
%   cnt_t: v_cnt*1 vector, store the ending index of each vertices
f_cnt = size(fout,1);
help_t = zeros(3*f_cnt,2);
help_t(1:f_cnt,1) = fout(:,1);
help_t(1:f_cnt,2) = [1:f_cnt]';
help_t(f_cnt+1:2*f_cnt,1) = fout(:,2);
help_t(f_cnt+1:2*f_cnt,2) = [1:f_cnt]';
help_t(2*f_cnt+1:3*f_cnt,1) = fout(:,3);
help_t(2*f_cnt+1:3*f_cnt,2) = [1:f_cnt]';

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
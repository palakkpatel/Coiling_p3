function grid_location = update_grid(p,v_min,v_max,grid_n)
% compute the grid location of point
% Input:
%   p: 1*3 vector. point
%   v_min: 1*3 vector
%   v_max: 1*3 vector
% Ouput:
%   grid_location: 1*3 vector
    grid_location = ceil((p-v_min)./((v_max-v_min)/grid_n));
end
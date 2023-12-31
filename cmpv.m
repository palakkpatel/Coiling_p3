function result = cmpv(x, y)
    %Comparing two vectors
% if ( size(x,2) == 3 && size(x,1) == 1 ...
%     || size(x,1) == 3 && size(x,2) == 1) ...
%     && (size(y,2) == 3 && size(y,1) == 1 ...
%     || size(y,1) == 3 && size(y,2) == 1)
    result = ( x(1) == y(1) & x(2) == y(2) & x(3) == y(3));
% else
%     error('x, y should be an array of size 1*3 or 3*1');
% end
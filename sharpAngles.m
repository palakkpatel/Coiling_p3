function [a_d] = sharpAngles(meshpoint, preshape)
% compute the bending energy of the coil by computing change in angles
% Input:
%   meshpoint: n*3 matrix, each row is a vertex of the coil.
%   preshape: n*3 matrix, each row is a vertex of the preshape.
%   sp: starting point
%   ep: ending point
% Output:
%   be: scaler, bending energy of the coil
%   g: n*3 matrix, return the gradient if needed

% meshpoint = [sp;meshpoint;ep];
len = size(meshpoint,1);
a_d=zeros(len,1);
sina = zeros(size(meshpoint,1),1); % sin(phi-phi')
cosa = zeros(size(meshpoint,1),1);
tana = zeros(size(meshpoint,1),1);
sinp = zeros(size(meshpoint,1),1); % sin(phi) of preshape
cosp = zeros(size(meshpoint,1),1); % cos(phi) of preshape

for ii = 2:len-1
    m1 = meshpoint(ii-1,:)-meshpoint(ii,:);
    m2 = meshpoint(ii+1,:)-meshpoint(ii,:);
    p1 = preshape(ii-1,:)-preshape(ii,:);
    p2 = preshape(ii+1,:)-preshape(ii,:);
    li = norm(m1) + norm(m2);
    cosp(ii) = p1*p2'/norm(p1)/norm(p2);
    a_p=acos(cosp(ii));
    cosm = m1*m2'/norm(m1)/norm(m2);
    a_m=acos(cosm);
    a_d(ii)=abs(a_p-a_m);
end

end
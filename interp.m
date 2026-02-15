% ---------------------------------------------------------------------- %
%   interp.m
%
%   Interpolates an Eulerian grid velocity field u onto Lagrangian marker
%   positions X. Periodicity is handled by wrapping grid indices with mod.
%
%   Inputs:
%       u - velocity field on grid, size [N x N x 2]
%       X - marker positions, size [n x 2]
%
%   Output:
%       U - interpolated velocities at markers
% ---------------------------------------------------------------------- %

function U = interp(u,X)
global h N;   

n = length(X(:,1));
U = zeros(n,2);

for k = 1:n          
  s = X(k,:)/h; % the real X moves in an infinitely large (i.e. not periodic) domain
  i = floor(s);           
  r = s-i;                
  i1 = mod((i(1)-1):(i(1)+2),N)+1;
  i2 = mod((i(2)-1):(i(2)+2),N)+1;
  w = phi1(r(1)).*phi2(r(2)); 
  U(k,1) = sum(sum(w.*u(i1,i2,1)));
  U(k,2) = sum(sum(w.*u(i1,i2,2))); 
end


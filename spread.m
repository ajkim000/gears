% ---------------------------------------------------------------------- %
%   spread.m
%
%   Spreads force densities from IB markers onto the Eulerian
%   fluid grid. 
%
%   Inputs
%       F - force densities at IB markers
%       X - massless IB marker positions
%
%   Output
%       force field on fluid grid
% ---------------------------------------------------------------------- %

function f = spread(F,X)
global h N dA Nb;

c = dA/(h*h);
f = zeros(N,N,2);

for k = 1:Nb
  s = X(k,:)/h;
  i = floor(s);
  r = s-i; 
  i1 = mod((i(1)-1):(i(1)+2),N)+1;        % 1 x 4
  i2 = mod((i(2)-1):(i(2)+2),N)+1;
  w = phi1(r(1)).*phi2(r(2));             % 4 x 4 
  f(i1,i2,1) = f(i1,i2,1)+(c*F(k,1))*w;   % 4 x 4 
  f(i1,i2,2) = f(i1,i2,2)+(c*F(k,2))*w;
end
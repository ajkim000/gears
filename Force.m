% Input
%   X: vector of IB points (size Nb x 2) 
% Output
%   F: vector of force densities (size Nb x 2)

%   Each pair of values refers to force felt in a small region (hence it
%   being a force DENSITY) near the corresponding IB point. 
%   Thus F*dA is the force transmitted to the fluid by volume element dA of the body.

function F = Force(X,Y)
global K;
F = K * (Y - X); % Restorative forces from springs
% ---------------------------------------------------------------------- %
%   Torque.m
%
%   Computes the torque on passive gear (gear 2) from force densities. 
%
%   Inputs
%       F2 - force densities from gear 2 markers
%       Y2 - massive boundary of gear 2
%       Y2_cm - center of mass of gear 2
% ---------------------------------------------------------------------- %

function T = Torque(F2,Y2,Y2_cm)
global dA 

T = 0;
for k = 1:length(Y2(:,1))
    T = T + dA * det([Y2(k,:) - Y2_cm; -F2(k,:)]);
end
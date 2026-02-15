% ---------------------------------------------------------------------- %
%   Force.m
%
%   Computes penalty spring forces enforcing rigidity in the pIB method.
%   Each Lagrangian marker is tethered to its corresponding rigid-body
%   position via a linear spring with stiffness K.
%
%   Inputs:
%       X - massless IB marker positions
%       Y - corresponding rigid-body positions
%
%   Output:
%       F - force density transmitted to the fluid at each marker
% ---------------------------------------------------------------------- %

function F = Force(X,Y)
global K;

F = K * (Y - X);
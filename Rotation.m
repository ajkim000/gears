% ---------------------------------------------------------------------- %
% Rotation.m
%
% Returns the 2x2 rotation matrix for a given angle theta (radians).
% ---------------------------------------------------------------------- %

function R = Rotation(theta)
  
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
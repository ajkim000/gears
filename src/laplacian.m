% ---------------------------------------------------------------------- %
% laplacian.m
%
% Computes the 2D finite-difference Laplacian of a vector field u
% on a periodic grid using a standard 5-point stencil.
% ---------------------------------------------------------------------- %

function w = laplacian(u)
global im ip h;

w = (u(ip,:,:)+u(im,:,:)+u(:,ip,:)+u(:,im,:)-4*u)/(h*h);


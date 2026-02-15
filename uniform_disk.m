% ----------------------------------------------------------------------%
% This function produces a set of points that form a uniform            %
% distribution on a disk of given radius.                               %
%                                                                       %
% Inputs                                                                %     
%   R – radius of the disk to be discretized                             %
%   N_rings – number of rings with which to discretize the disk          %
%                                                                       %
% Notes                                                                 %
%   Radial spacing will be slightly different from arclength spacing    %
% ----------------------------------------------------------------------% 


function [x, y, dA] = uniform_disk(R, N_rings)

    count = 0;
    x = zeros(1, 3*N_rings^2);
    y = zeros(1, 3*N_rings^2);
    
    for n = 1:N_rings
        r = (n-0.5) * (R/N_rings);      % radius of ring n
        m = 3 * (2*n - 1);              % number of points on ring n
        dtheta = 2*pi/m;
        theta = pi * mod(n,2);          % start at pi if n is odd, 0 if n is even
        for mm = 1:m
            count = count + 1;
            x(count) = r * cos(theta);
            y(count) = r * sin(theta);
            theta = theta + dtheta;
        end
    end

    x = x.';
    y = y.';
    dA = pi*R^2 / count;
end

% ---------------------------------------------------------------------- %
%   initialize.m
%
%   Sets up parameters, geometry, and initial state for the 2D
%   penalty immersed boundary simulation of two rotating gears.
%
%   Initializes:
%       - Physical and numerical parameters
%       - Massless and massive IB boundaries for gears and container
%       - Rigid-body frames for gear rotation
%       - Fluid velocity field
%       - Particle visualization
%
%   Units: cm, g, s.
% ---------------------------------------------------------------------- %

% -------------------------------------------- %
%             INITIAL CONSTANTS                %
% -------------------------------------------- %

% USER PARAMETERS (edit these)

target = '../trials/trial_0';   % where png frames and mat files are saved 
Omega1 = -2*pi;              % active gear angular velocity (rad/s)
gap    = 1.25*2.54;          % edge-to-edge interior gap between gears (cm)
K      = 100;                % penalty stiffness
N      = 256;                % Eulerian grid resolution (NxN)
dt     = 0.01;               % time step (s)
tmax   = 10;                 % final simulation time (s)
clockmax = ceil(tmax/dt);    

% Output cadence
num_img = 500;               % number of png frames
num_mat = 5;                 % number of mat snapshots

% Particle viz
Np = 1500;                 
marker_size = 10;
tail_len = 1000;     

% Create output directory
if exist(target,'dir'); rmdir(target,'s'); end
mkdir(target)

% Physical parameters
rho = 1.2012;               % fluid density (g/cm^3)
mu = 2.4024;                % fluid dynamic viscosity (g/(cm s))
R = 1.905;                  % gear radius (cm)
D = 17.145;                 % container diameter (cm)
domain = 1.2*D;             % domain should be large enough to fit the container

% Grid
h = domain/N;               % grid spacing
ip = [(2:N),1];         
im = [N,(1:(N-1))];  
xgrid = zeros(N,N);
ygrid = zeros(N,N);
for j=0:(N-1)
    xgrid(j+1,:)=j*h;
    ygrid(:,j+1)=j*h;
end
mesh_ratio = 1;             % {solid mesh size} / {fluid mesh size} 

% Initialize state of passive gear
Omega2 = 0;                 % initial angular velocity on gear 2 (rad/s)
Omega2_all = [];
T = 0;                      % initial torque on gear 2
Torque_all = [];
L = 0;                      % initial angular momentum on gear 2
I0 = 0.7450;                % moment of inertia of gear 2 (g cm^2)
num_layers = 10;            % number of layers in the container to prevent leakage

% -------------------------------------------- %
%         INITIALIZE MASSLESS BOUNDARY         %
% -------------------------------------------- %

% Initialize IB points on gears
h_gear_desired = mesh_ratio * h;                % desired IB spacing on the gears (ideally similar to h)
N_rings = ceil(R/h_gear_desired);               % radial spacing ~ R/N (see uniform_disk.m)
h_gear = R/N_rings;                             % note h_gear will not be exactly equal to h_gear_desired
[x,y,dA] = uniform_disk(R, N_rings);            % IB points on gears
X = [];
X(:,1) = [x + (domain/2) - R - gap/2; x + (domain/2) + R + gap/2];    
X(:,2) = [y + domain/2; y + domain/2];                                  

% Initialize IB points on container
N_rings_temp = ceil((D/2) / h_gear);            % number of rings needed to achieve a spacing similar to h_gear
h_cont = (D/2) / N_rings_temp;                  % actual radial spacing of the container
R_cont_outer = (D/2) + num_layers * h_cont;     % outer radius of the container can be determined from h_cont
[z1,z2,dA_cont] = uniform_disk(R_cont_outer,N_rings_temp + num_layers); % check: dA_cont should be very similar to dA
dist = sqrt(z1.^2 + z2.^2);                     % distance of each point from origin
mask = dist >= (D/2);                           % mask to keep only those in the annulus
z1 = z1(mask);
z2 = z2(mask);
Z = domain/2 + [z1,z2];                         % translate container center to domain center

% Concatenate gear and container boundaries
Nb_gears = length((X(:,1)));                    % total number of IB points on both gears
Nb_cont = length((Z(:,1)));                     % IB points on container
X = [X;Z];                                      % concatenate boundaries
Nb = Nb_gears + Nb_cont;                        % total IB points

% -------------------------------------------- %
%        INITIALIZE MASSIVE BOUNDARY           %
% -------------------------------------------- %

Y = X;

% Gear 1
Y1 = Y(1:Nb_gears/2,:);                     % massive boundary
Y1_cm = [domain/2 - R - gap/2, domain/2];   % center of mass
s1_end = (Y1_cm + [R, 0]).';                % radius marker to track rotation

% Gear 2
Y2 = Y(Nb_gears/2+1:Nb_gears,:);            % massive boundary
Y2_cm = [domain/2 + R + gap/2, domain/2];   % center of mass
s2_end = (Y2_cm + [R, 0]).';                % radius marker to track location

% -------------------------------------------- %
%        INITIALIZE ORTHONORMAL FRAME          %
% -------------------------------------------- %

% Gear 1
C1 = Y1 - Y1_cm;                            % relative coordinates (wrt center of mass)
E1 = eye(2,2);                              % initial orthonormal frame

% Gear 2
C2 = Y2 - Y2_cm;                            % relative coordinates (wrt center of mass)
E2 = eye(2,2);                              % initial orthonormal frame

% -------------------------------------------- %
%          INITIALIZE FLOW FIELD               %
% -------------------------------------------- %

u = zeros(N,N,2);                           % initially static fluid flow

% -------------------------------------------- %
%                 ANIMATION                    %
% -------------------------------------------- %

%% Particle visualization

% Seed domain with particles
Xp = domain * rand(Np, 2);                  % Np randomly generated particles
colors = zeros(Np, 3); above = Xp(:,2) > domain/2; below = ~above;
color_scheme = orderedcolors("glow12"); 
colors(above,:) = repmat(color_scheme(1,:), sum(above), 1);
colors(below,:) = repmat(color_scheme(2,:), sum(below), 1);

% Only keep particles outside gears and inside container
dist_gear1 = sqrt((Xp(:,1) - Y1_cm(1)).^2 + (Xp(:,2) - Y1_cm(2)).^2);
dist_gear2 = sqrt((Xp(:,1) - Y2_cm(1)).^2 + (Xp(:,2) - Y2_cm(2)).^2);
dist_center = sqrt((Xp(:,1) - (domain/2)).^2 + (Xp(:,2) - (domain/2)).^2);
keep = (dist_gear1 > R) & (dist_gear2 > R) & (dist_center < D/2);

% Initialize particle history matrix
Xp_all = nan(Np, 2, tail_len+1);
Xp_all(:,:,1) = Xp; 

% Plot massive boundary and particles for sanity check
figure;
set(gcf,'double','on')
hold on
plot(Y(:,1),Y(:,2),'ko', 'MarkerSize', 5) % plot massive boundary
plot([Y1_cm(1), s1_end(1)],[Y1_cm(2), s1_end(2)],'g-','LineWidth',3)
plot([Y2_cm(1), s2_end(1)],[Y2_cm(2), s2_end(2)],'g-','LineWidth',3)
scatter(Xp(keep,1), Xp(keep,2), marker_size, colors(keep,:), 'filled');
axis equal
axis manual
axis off
axis([0,domain,0,domain])
set(gca,'fontsize',16)
drawnow
hold off
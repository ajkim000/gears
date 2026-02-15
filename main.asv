% ------------------------------------------------------------------------- %
%   Hydrodynamic Spin-Coupling of Two Rotors (2D)
%
%   This script simulates the fluid-mediated coupling between two rigid,
%   disk-shaped rotors immersed in a viscous incompressible fluid using the
%   rigid penalty immersed boundary (pIB) method.
%
%   The left rotor ('active') is driven at a prescribed angular velocity.
%   The emergent angular velocity of the right rotor ('passive') is 
%   measured to quantify fluid-mediated spin coupling.
%
%   Outputs
%       - Time-resolved passive angular velocity (Omega2_all)
%       - Hydrodynamic torque history (Torque_all)
%       - Gear and particle visualization frames (png files)
%       - Periodic simulation snapshots (.mat files)
%
%   Requirements
%       Helper functions must be on the MATLAB path:
%       fluid.m, Force.m, init_a.m, initialize.m, interp.m, laplacian.m,
%       phi1.m, phi2.m, Rotation.m, sk.m, skew.m, spread.m, Torque.m,
%       uniform_disk.m 
%
%   Author: Alison Kim
%   Affiliation: Courant Institute of Mathematical Sciences, NYU
%   Date: February 2026
%
%   References
%       Kim, Y. & Peskin, C.S. (2016).
%           "A penalty immersed boundary method for a rigid body in fluid."
%           Physics of Fluids, 28(3), 033603.
%
%       Battista, N.A., Strickland, W.C., & Miller, L.A. (2017).
%           "IB2d: a Python and MATLAB implementation of the immersed 
%           boundary method." Bioinspiration & Biomimetics, 12(3), 036003.
% ------------------------------------------------------------------------- %

global dt Nb Nb_gears N h rho mu ip im a K dA
initialize 
init_a 
img_interval = floor(clockmax / num_img);   % how often to save frame images
mat_interval = floor(clockmax / num_mat);   % how often to save mat file
fig = figure;

%% Initial animation

clock = 0;

% Parametrize boundary
th = linspace(0, 2*pi, 200);
gear1 = [-(R+gap/2) + domain/2 + R*cos(th); (domain/2) + R*sin(th)];
gear2 = [R+gap/2 + domain/2 + R*cos(th); (domain/2) + R*sin(th)];
cont_outer = (domain/2) + (1.05*D/2) .* [cos(th); sin(th)];
cont_inner = (domain/2) + (D/2) .* [cos(th); sin(th)];

% Draw boundary and particles
clf(fig)
hold on
plot(gear1(1,:), gear1(2,:), 'k-', 'LineWidth', 2)
plot(gear2(1,:), gear2(2,:), 'k-', 'LineWidth', 2)
plot(cont_outer(1,:), cont_outer(2,:), 'k-', 'LineWidth', 2)
plot(cont_inner(1,:), cont_inner(2,:), 'k-', 'LineWidth', 2)
plot([Y1_cm(1), s1_end(1)],[Y1_cm(2), s1_end(2)],'k-','LineWidth',2)
plot([Y2_cm(1), s2_end(1)],[Y2_cm(2), s2_end(2)],'k-','LineWidth',2)
scatter(Xp(keep,1), Xp(keep,2), marker_size, colors(keep,:), 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0); 
text(0.5, -0.08, sprintf('time = %.2f s', clock*dt), ...
    'Units','normalized', 'HorizontalAlignment','center', 'FontSize',14);
axis equal
axis off
drawnow
filename1 = sprintf('%s/f-%05.2f.png',target,clock*dt); % time label 
print(filename1,'-dpng','-r250') 
hold off

while clock < clockmax

  clock = clock + 1;
  
  % -------------------------------------------- %
  %   PRELIMINARY SUBSTEP (step forward dt/2)    %
  % -------------------------------------------- %
  
  %% Update massless boundary
  XX = X + (dt/2)*interp(u,X);     
  XXp = Xp + (dt/2)*interp(u,Xp);       % update particles 

  %% Update massive boundary

  % Gear 1
  theta1_half = Omega1 * (dt/2);
  EE1 = Rotation(theta1_half) * E1;     % update orthogonal frame according to Omega1
  YY1 = Y1_cm + (EE1*C1.').';           % update massive boundary (rotate)

  % Gear 2
  Omega2 = L/I0;                        % new angular velocity of gear 2
  theta2_half = Omega2 * (dt/2);
  EE2 = Rotation(theta2_half) * E2;     % new orthonormal frame
  YY2 = Y2_cm + (EE2*C2.').';           % update massive boundary
  YY = [YY1; YY2; Z];
  
  %% Spread forces and move fluid
  FF = Force(XX,YY);                    % F[dA] = force by solid element dA on nearby fluid
  ff = spread(FF,XX);                   % f[x,y] = force by whole body on fluid element [x,y]
  [u,uu] = fluid(u,ff);                 % update fluid flow due to spreaded force
   FF2 = FF(Nb_gears/2+1:Nb_gears,:);
  TT = Torque(FF2,YY2,Y2_cm);           % compute torque on gear 2
  Torque_all = [Torque_all TT];
  LL = L + (dt/2) * TT;                 % update angular momentum of gear 2
 
  % -------------------------------------------- %
  %        FINAL SUBSTEP (step forward dt)       %
  % -------------------------------------------- %

  % Update massless boundary
  X = X + dt*interp(uu,XX); 

  % Update particles and their tails
  Xp = Xp + dt*interp(uu,XXp);              % update particle positions
  Xp_all(:,:,2:end) = Xp_all(:,:,1:end-1);  % bump off oldest particles
  Xp_all(:,:,1) = Xp;                       % push new particles into tail

  % Update orthonormal frames
  theta1 = Omega1 * dt;
  E1 = Rotation(theta1) * E1;       % update gear 1 orthonormal frame
  Omega2 = LL/I0;                   % gear 2 new angular velocity
  Omega2_all = [Omega2_all Omega2];
  theta2 = Omega2 * dt;
  E2 = Rotation(theta2) * E2;       % update gear 2 orthonormal frame
   
  % Update massive boundary
  Y1 = Y1_cm + (E1*C1.').';         % gear 1
  Y2 = Y2_cm + (E2*C2.').';         % gear 2
  Y = [Y1;Y2;Z];

  % Update angular momentum for next step
  L = L + dt*TT;

  % -------------------------------------------- %
  %                   ANIMATION                  %
  % -------------------------------------------- %

  % Update gear radii markers to show rotation
  s1_end = Rotation(theta1) * (s1_end - Y1_cm.') + Y1_cm'; % gear 1
  s2_end = Rotation(theta2) * (s2_end - Y2_cm.') + Y2_cm'; % gear 2

  % Plot frame
  if mod(clock, img_interval) == 0

      clf(fig); 

      % Draw boundary
      plot(gear1(1,:), gear1(2,:), 'k-', 'LineWidth', 2)
      hold on
      plot(gear2(1,:), gear2(2,:), 'k-', 'LineWidth', 2)
      plot(cont_outer(1,:), cont_outer(2,:), 'k-', 'LineWidth', 2)
      plot(cont_inner(1,:), cont_inner(2,:), 'k-', 'LineWidth', 2)
      plot([Y1_cm(1), s1_end(1)],[Y1_cm(2), s1_end(2)],'k-','LineWidth',2)
      plot([Y2_cm(1), s2_end(1)],[Y2_cm(2), s2_end(2)],'k-','LineWidth',2)
  
      % Draw particles and their tails 
      scatter(Xp(keep,1), Xp(keep,2), marker_size, colors(keep,:), 'filled', ...
          'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0); 
      alphas = linspace(0.2,0.01,tail_len);         % tail transparency
      for tail_idx = 1:tail_len  
          pts = Xp_all(keep,:,tail_idx+1);     
          tail_marker_size = (-marker_size/tail_len) * (tail_idx-1) + marker_size;
          scatter(pts(:,1), pts(:,2), tail_marker_size, colors(keep,:), 'filled',...
              'MarkerFaceAlpha', alphas(tail_idx), 'MarkerEdgeAlpha', 0);
      end
      axis equal
      axis off

      % Time label
      text(0.5, -0.08, sprintf('time = %.2f s', clock*dt), ...
         'Units','normalized', 'HorizontalAlignment','center','FontSize',14);
      hold off

      drawnow
      filename1 = sprintf('%s/f-%05.2f.png',target,clock*dt); 
      print(filename1,'-dpng','-r250') 
      hold off
  end

  % Save data
  if mod(clock, mat_interval) == 0
      filename2 = sprintf('%s/s-%05.2f.mat',target,clock*dt);
      save(filename2,'u','L','Omega1','Omega2','Torque_all','Omega2_all', ...
          'X','Y','Y1_cm','Y2_cm','s1_end','s2_end','Z','C1','C2','E1','E2',...
          'gap','K','N','dt','I0','mesh_ratio', ...
          'rho','mu','R','D','domain','h','ip','im','Nb','Nb_gears','dA', ...
          'xgrid','ygrid','clock',...
          'Xp','Xp_all','colors','keep','marker_size','tail_len') % particle info
  end
  hold off

end
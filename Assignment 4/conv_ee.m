function [phi_out, phi_a_out, err_mean, CFL]=conv_ee(U0,xend,points,tsteps,dt)

% Variable allocation
dx = xend / ( points - 1 );
x = 0.0 : dx : xend;
tend   = dt * tsteps;
CFL = (U0*dt)/dx;
phi_out(1,:) = x; % x values for plotting
phi_out_a(1,:) = x; % x values for plotting
k = 2;

% Initial condition
 phi = sin(x);

% % Check initial field:
% plot(x, phi, 'r');
% title("initial field")
    
%% Explicit Euler

    for i = 1 : tsteps

      % Periodic boundary conditions at x=0:
      phinew(1) = phi(points-1);

      % Loop over grid points in space:
      for j = 2 : points - 1

        phinew(j) = phi(j) - U0*((phi(j+1)-phi(j-1))/(2*dx))*dt;

      end

      % Periodic boundary conditions at x=2*pi:
      phinew(points) = phi(2);

      % Write new field back to old field:
      phi = phinew;

      % Analytical solution
%       for j = 1 : points
%           phi_a(j)=sin(x(j)-U0*(i*dt));
%       end
    phi_a = sin(x-U0*(i*dt));
%     phi_a = exp(x-U0*(i*dt));
      
%       if ~mod(i,10) == 1
%           % Plot transported wave for each timestep
%           figure(2)
%           plot(x, phi, 'r', x, phi_a, 'k');          
%           pause(0.003);
%       end
      title("animation")

      if i == 1 || i == tsteps/10 || i == tsteps/2 || i == tsteps
          phi_out(k,:) = phi;
          phi_a_out(k,:) = phi_a;
          k = k+1;
      end
      err_mean(i,1) = i;
      err_mean(i,2) = sqrt(mean((phi_a+999 - phi+999).^2))/mean(phi_a+999);
    end

end 



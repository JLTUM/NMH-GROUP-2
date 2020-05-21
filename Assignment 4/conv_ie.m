function [phi_out, phi_a_out, err_mean, CFL]=conv_ee(U0,xend,points,tsteps,dt)

% Variable allocation
dx = xend / ( points - 1 );
x = 0.0 : dx : xend;
tend   = dt * tsteps;
CFL = (U0*dt)/dx;
phi_out(1,:) = x; % x values for plotting
phi_out_a(1,:) = x; % x values for plotting
k = 2;
A   = zeros(points,points);
b   = zeros(points,1);
phi = zeros(points,1);

% Initial Condition
phi = sin(x)';

% Check initial field:
plot(x, phi, 'r');
hold on;

%% Implicit Euler

% Compute coefficients of matrix A
a_w = -(U0*dt)/(2*dx);
a_p = 1;
a_e = (U0*dt)/(2*dx);

for i = 1 : tsteps

 % Periodic boundary conditions at x=0:
  A(1,1) = a_p;
  A(1,2) = a_e;
  A(1,points-1) = a_w;
  b = phi;

  % Loop over grid points in space:
  for j = 2 : points - 1
      
      A(j,j) = a_p;
      A(j,j+1) = a_e;
      A(j,j-1) = a_w;

  end

  % Periodic boundary conditions at x=2*pi:
  A(points,points) = a_p;
  A(points,2) = a_e;
  A(points,points-1) = a_w;
    
  % Solve the linear system of equations
  phi = (A\b);
  
  % Analytical solution
%   for j = 1 : points
%           phi_a(j)=sin(x(j)-U0*(i*dt));
%   end
  
  phi_a = sin(x - U0*(i*dt));

    if i == 1 || i == tsteps/10 || i == tsteps/2 || i == tsteps
      phi_out(k,:) = phi;
      phi_a_out(k,:) = phi_a;
      k = k+1;
    end
      err_mean(i,1) = i;
      err_mean(i,2) = sqrt(mean((phi_a+999 - phi+999).^2))/mean(phi_a+999);

  % Plot transported wave for each timestep
  plot(x, phi, 'r', x, phi_a, 'k');
  hold off;
  pause(0.003);

end
toc

disp(toc)

end



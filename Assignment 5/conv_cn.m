function [phi_out, phi_a_out, CFL, Pe_cell]=conv_cn(Gamma,U0,points,dt)

% Variable allocation
xend   = 2.0 * pi;
dx = xend / ( points - 1 );
x = 0.0 : dx : xend;
tsteps = 1000;
CFL = (U0*dt)/dx;
Pe_cell = (U0*dx)/Gamma;
phi_out(1,:) = x; % x values for plotting
phi_a_out(1,:) = x; % x values for plotting

% Preallocation
A       = zeros(points,points);
A_bar   = zeros(points,points);

% Initial Solution
phi = sin(x)';

%% Compute coefficients of matrix A

    a_w = U0/(2*dx)+Gamma/(dx^2);
    a_p = -(2*Gamma)/(dx^2);
    a_e = -U0/(2*dx)+Gamma/(dx^2);
    
    % Coefficients of unknown timestep
    A_w = (-dt/2*a_w);
    A_p = 1-(dt/2)*a_p;
    A_e = (-dt/2*a_e);
    
    % Coefficients of known timestep
    A_w_bar = -A_w;
    A_p_bar = 1+(dt/2)*a_p;
    A_e_bar = -A_e;

%% Allocation of A and A_bar

    % Periodic boundary conditions at x=0:
    A(1,1) = A_p;
    A(1,2) = A_e;
    A(1,points-1) = A_w;
    
  % Loop over grid points in space:
  for j = 2 : points - 1
    A(j,j) = A_p;
    A(j,j+1) = A_e;
    A(j,j-1) = A_w;
  end

  % Periodic boundary conditions at x=2*pi:
    A(points,points) = A_p;
    A(points,2) = A_e;
    A(points,points-1) = A_w;
    
  % A-Bar allocation
    A_bar(A==A_p) = A_p_bar;
    A_bar(A==A_e) = A_e_bar;
    A_bar(A==A_w) = A_w_bar;
       
%% Loop over timesteps
for i = 1 : tsteps
    
  % Assign matrix of known timestep
    b = A_bar*phi;

  % Solve the linear system of equations
    phi = A\b;

  % Analytical solution 
    phi_a = exp(-Gamma*(i)*dt) * sin(x - U0*((i)*dt));

  % write solution to matrix for plot
    if ~mod(i,2) == 1
        phi_out(end+1,:) = phi;
        phi_a_out(end+1,:) = phi_a;
    end
    
%   % Plot transported wave for each timestep
%     plot(x, phi, 'r', x, phi_a, 'g');
%     hold off;
%     pause(0.03);

end

end


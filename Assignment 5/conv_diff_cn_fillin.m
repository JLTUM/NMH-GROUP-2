function [phi_out, phi_a_out, CFL]=conv_ee(U0,points,dt)

close all
format long;
clear;
hold off;

% Set convection velocity
U0 = 1;

% Set diffusion coefficient
Gamma = 1.0;

% Discrete spacing in space
xend   = 2.0 * pi;
points = 40; 
dx     = xend / ( points - 1 );
% Grid with x locations:
x = 0.0 : dx : xend;

% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
dt     = 0.1;
tend   = dt * tsteps;

% Initialise coefficient matrix A, constant vector b
% and solution vector phi
A   = zeros(points,points);
A_bar   = zeros(points,points);
b   = zeros(points,1);
phi = zeros(points,1);

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

  % Plot transported wave for each timestep
    plot(x, phi, 'r', x, phi_a, 'g');
    hold off;
    pause(0.03);

end

end


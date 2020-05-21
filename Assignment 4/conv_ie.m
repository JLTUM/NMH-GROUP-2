% Course: CFD Lab
% TU Muenchen, Summer term 2020
%
%Group 2
%Assignment 4 - Main code
%Andreas Mirlach, Julian Lenz, Faro Schäfer, Nick Pfeiffer
%
% This is the Matlab script for the unsteady 1D convection equation
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi
%  ------- =  -U0 * -------
%    dt               dx
%
% 0 <= x <= 2pi
%
% periodic boundary condition
% phi(0) = phi(2pi)
%
% initial condition
% t = 0  ==>  phi = sin(x)
%
% Central difference scheme (CDS) for spatial discretization
% Implicit Euler scheme for time advancement
%
%--------------------------------------------------------------------------
%
clear ,clc, close all;
% function [matrix]=conv_ee(U0,xend,points,tsteps,
% Set convection velocity
U0 = 1.0;

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
b   = zeros(points,1);
phi = zeros(points,1);

% Initialise the solution (initial condition)
% Loop over grid points in space:
for j = 1 : points
   % phi(j) = some function of x(j) ???;
end

% Check initial field:
plot(x, phi, 'r');
hold on;
pause(3);

% Compute coefficients of matrix A
% a_w = ???
% a_p = ???
% a_e = ???

% Implicit Euler:
%----------------
%
% Loop over timesteps:
for i = 1 : tsteps

  % Periodic boundary conditions at x=0:
  % A(???,???) = ???
  % A(???,???) = ???
  % A(???,???) = ???
  % b(???) = ???

  % Loop over grid points in space:
  for j = 2 : points - 1

    % A(???,???) = ???
    % A(???,???) = ???
    % A(???,???) = ???
    % b(???) = ???

  end

  % Periodic boundary conditions at x=2*pi:
  % A(???,???) = ???
  % A(???,???) = ???
  % A(???,???) = ???
  % b(???) = ???

  % Solve the linear system of equations
  phi = A\b;

  % Analytical solution
  for j = 1 : points
    % phi_a(j) = some function of x(j) & t(i) ???
    % hint: t(i) = i * dt
  end

  % Plot transported wave for each timestep
  plot(x, phi, 'r', x, phi_a, 'g');
  hold off;
  pause(0.003);

end


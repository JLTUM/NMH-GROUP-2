% Course: NHM-Lab
% SS 2020

% This is the Matlab script for the convection-diffusion equation 

% You must fill in the missing parts by yourself!
% Missing parts are marked by ???

% Tianshi Sun
% tianshi.sun@tum.de


% solution of

%       d phi          d^2 phi
%  -U0 ------- + Gamma ------- = 0
%        dx             dx^2

% u=0
% 0 <= x <= 2pi

% Clear all variables and plots.
format long;
clear;
hold off;

% Set convection velocity
U0 = 1.0;

% Set diffusivity
Gamma = 1.0;

% Discrete spacing in space!
xend   = 2.0*pi;
points = 40; 
dx     = xend/(points-1);
% Grid with x locations:
x = 0.0 : dx : xend;

% Initialization of field
phi = zeros(points,1);

% Initialization of matrix A
A = zeros(points,points);

% Initialization of vector b
b = zeros(points,1);

% Boundary condition
phi_0   = 0.0;
phi_end = 1.0;

% Loop over grid points in space
% note that boundary points are excluded
for i = 2 : points-1

     a_w = ???
     a_p = ???
     a_e = ???
     
%     assign values to matrix A


%     assign values to vector b



end

% Boundary conditions

% at i = 1
 A(???,???) = ???
 b(???) = ???

% at i = points
 A(???,???) = ???
 b(???) = ???

% Solution of the linear system
phi = A\b

%Analytical solution

phi_analytic = ???

%error

nn = points / 2;
er = ???

% Plot the solution
plot(x,phi,'r', x,phi_analytic, 'go');

% Plot the error

???
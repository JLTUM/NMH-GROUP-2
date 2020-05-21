%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 1
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

% ----------------------------------------------------------------------- %

% This is a Matlab code to solve 1D steady advection-diffusion equation
% discritised by finite-volume schemes 
%
% differential form of advection-diffusion eqn.:
%
%       d phi          d^2 phi
%  -U0 ------- + Gamma ------- = 0
%        dx             dx^2
%
% You are asked to fill in the missing parts to complete the implementation.
% Missing parts are marked by ???
%
% Author: Yoshiyuki Sakai
% Email: yoshiyuki.sakai@tum.de

function A_D_eq(U0,Gamma,points,nn,scheme)
% Clear all variables and plots.
format long;
clear;
hold off;

% Set constant advection velocity
U0 = 1.0;

% Set constant diffusivity
Gamma = 1.0;

% Set up grid cells
xend = 2.0 * pi;
cells = 51; 
dx = ???

% Array of grid cell centre locations:
x = ??? 

% Initialization of cell-averaged field
phi = zeros(cells,1);

% Initialization of matrix A
A = zeros(cells,cells);

% Initialization of vector b
b = zeros(cells,1);

% Boundary cell face values
phi_0   = 0.0;
phi_end = 1.0;

% Loop over grid cells
% NB: boundary cells are excluded
for i = 2 : cells-1

     a_w = ???
     a_p = ???
     a_e = ???
     
%     assign values to LHS matrix A
     A(i,i) = a_p;
     A(i,i-1) = a_w;
     A(i,i+1) = a_e;

%     assign values to RHS vector b
     b(i) = ???
end

% Boundary conditions (Dirichlet at boundary cell faces)

% at i = 1
 A(1,1) = ???
 A(1,2) = ???
 b(1) = ???

% at i = cells
 A(cells,cells) = ???
 A(cells,cells-1) = ???
 b(cells) = ???

% Solution of the linear system
phi = A\b;

% Compute analytical solution
phi_analytic = ???

% Compute relative error
nn = ceil(cells / 2);
err_rel = ???

% Compute mean error
err_mean = ???

% Plot the numerical and analytical solutions
plot(x, phi, 'r', x, phi_analytic, 'go');

% Plot the error as function of dx in log-log scale
???

end
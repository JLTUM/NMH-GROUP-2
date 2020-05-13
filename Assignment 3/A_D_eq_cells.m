%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 3
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
%%

% Set constant advection velocity
U0 = 1.0;

% Set constant diffusivity
Gamma = 1.0;

% Set up grid cells
xend = 2.0 * pi;
cells = 51;

%%
%function A_D_eq_cells(U0,Gamma,cells,nn,scheme)

% Clear all variables and plots.
% format long;
% clear;
% hold off;

% Set up grid cells
xend = 2.0 * pi;
dx = xend/(cells-1);% points = 10;

% Array of grid cell centre locations:
x = 0.0 : dx : xend;

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

     a_w = U0/2+Gamma/dx;
     a_p = -2*Gamma/dx;
     a_e = -U0/2+Gamma/dx;
     
%     assign values to LHS matrix A
     A(i,i) = a_p;
     A(i,i-1) = a_w;
     A(i,i+1) = a_e;

%     assign values to RHS vector b
     b(i) = 0; %???
end

% Boundary conditions (Dirichlet at boundary cell faces)

% at i = 1
 A(1,1) = -U0/2-3*Gamma/dx ; %Ap,BC left side
 A(1,2) = -U0/2+Gamma/dx ; %Ae
 b(1) = (-U0-2*Gamma/dx)*phi_0; 

% at i = cells
 A(cells,cells) = U0/2-3*Gamma/dx ; %Ap,BC right side
 A(cells,cells-1) = U0/2+Gamma/dx ; %Aw
 b(cells) = U0-2*Gamma/dx*phi_end;

% Solution of the linear system
phi = A\b;

% Compute analytical solution
phi_analytic = (exp((U0.*x/Gamma))-1)/(exp((2*pi*U0)/Gamma)-1); %%correct?

% % Compute relative error
% nn = ceil(cells / 2);
% err_rel = ???
% 
% % Compute mean error
% err_mean = ???

% Plot the numerical and analytical solutions
plot(x, phi, 'r', x, phi_analytic, 'go');
legend('Numerical','Analytic')

% Plot the error as function of dx in log-log scale
%???

%end
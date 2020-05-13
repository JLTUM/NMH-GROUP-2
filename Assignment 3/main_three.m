clear,clc, close all;
% Set constant advection velocity
U0 = 1.0;

% Set constant diffusivity
Gamma = 1.0;

% Set up grid cells
xend = 2.0 * pi;
cells = 51; 

nn = cells/2;
scheme = "Central"; % "Central" or "Upwind"

A_D_eq_cells(U0,Gamma,cells,nn,scheme)
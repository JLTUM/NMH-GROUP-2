%**************************************************************************
% NMH Lab Summer Semester 2020, Assignment 6
%
% This code solves the following vectorial 2D advection problem 
%
% "TWOD_CONDIFF_VAR"
% du/dt = -(u0*du/dx + v0*du/dy) + nu*(d^2 u/dx^2 + d^2 u/dy^2) 
% dv/dt = -(u0*dv/dx + v0*dv/dy) + nu*(d^2 v/dx^2 + d^2 v/dy^2)
% where u_tr and v_tr are functions of space
%
% For spatial discretisation as central difference scheme is used
% and for time integration the explicit Euler time integration and a 
% Runge-Kutta method are available.
%
% Periodic boundary condition is applied both in x- and y-direction
%
% author: H. Zeng & L. Unglehrt & D.Quosdorf & Y.Sakai
% June, 2020
%**************************************************************************

clear all
close all

%% Initialize simulation
% read infile 
infilename = 'infile_2D_condiff_var.mat';
fprintf('infilename is: %s\n', infilename)

% build structures 'grid' and 'flow'
[grid, flow] = build_structs;
%[grid_e, flow_e] = build_structs;
fprintf('struct built\n')

% fill some fields of 'grid' and 'flow' with data from infile
[grid, flow] = set_params(grid, flow, infilename);
%[grid_e, flow_e] = set_params(grid, flow, infilename);
fprintf('parameters set\n')

% ---- Generate an equidistant grid -----------------------------------
[grid] = generate_grid(grid);
%[grid_e] = generate_grid(grid);
    
% initialisation of flow field
[flow] = set_initial_condition(grid, flow);
fprintf('flow field initialised\n')

%% Time integration
%%1 = euler 0 = Runge
[flow,flow_e,error] = time_step_rk(grid, flow);


hold off
figure(2)
plot(error)


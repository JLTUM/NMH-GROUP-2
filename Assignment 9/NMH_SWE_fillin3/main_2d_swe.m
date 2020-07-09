%**************************************************************************
% NMH Lab Summer Semester 2020, Assignment 9
%
% This code solves the 2D shallow water equations
%
% author: H. Zeng & L. Unglehrt
% July, 2020
%**************************************************************************
clear;
close all

%% Initialize simulation
% read infile 
% infilename = ???
fprintf('infilename is: %s\n', infilename)

% build structures 
[grid, run, constants, flow, bconds] = build_structs;
fprintf('struct built\n')

% fill some fields of 'grid' and 'flow' with data from infile
[grid, run, constants] = set_params(infilename);
fprintf('parameters set\n')

% Generate an equidistant grid 
[grid] = generate_grid(grid);    

% Set initial conditions 
run.t = 0;
[ flow ] = set_initial_condition( grid, flow );

%% ---- Create boundary conditions -----------------------------------------
[ bconds ] = set_boundary_conditions();

% ---- Setup of for time integration --------------------------------------
% Frequency of diagnostic output
% itdiag = ???

%% Time integration

fprintf('start time integration\n')
for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );

    % Diagonostic output
    if mod(itstep, itdiag) == 0
        % CFL number
        fprintf('%d : CFL number:         %e\n', itstep, ...
            compute_CFL_number(constants, grid, run.dt, flow.h, flow.hu, flow.hv));
        
%% Plot results
% TODO TODO TODO TODO TODO TODO TODO
    end
end

%**************************************************************************
% NMH Lab Summer Semester 2020, Assignment 7
%
% This code solves the 2D Shallow Water equations
%
% author: H. Zeng & L. Unglehrt
% June, 2020
%**************************************************************************
clear;
close all

%% Initialize simulation
% read infile 
infilename = 'infile_2D_swe_test.mat';
fprintf('infilename is: %s\n', infilename)

% build structures 
[grid, run, constants, flow, bconds] = build_structs;
fprintf('struct built\n')

% fill some fields with data from input file
[grid, run, constants] = set_params(infilename);
fprintf('parameters set\n')

% Generate an equidistant grid 
[grid] = generate_grid(grid);    

% Set initial conditions 
run.t = 0;
[ flow ] = set_initial_condition( grid, flow );

% Create boundary conditions
bconds.bwest = {'WALL'};
bconds.beast = {'WALL'};
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

%% Time integration
for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );
    
    
% if itstep == 1
%     
    

%% Plot results
surf(grid.x, grid.y, flow.h+flow.zb,'FaceAlpha',0.5)
hold on 
surf(grid.x, grid.y, flow.zb,'FaceColor','b')
zlim([-1 2])
title(itstep)
pause(0.11)
hold off



end


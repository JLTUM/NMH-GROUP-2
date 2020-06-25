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
fprintf('grid set\n')

% Set initial conditions 
run.t = 0;
[ flow ] = set_initial_condition( grid, flow );

% Create boundary conditions
bconds.bwest = {'WALL'};
bconds.beast = {'WALL'};
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

% Compute and print diffusion number
 %   D = k_f*dt/(S_0*min(dx,dy)^2);
%fprintf('Diffusion number: %d\n',D);

%% Time integration
for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );
    
%% Plot results
%  figure(1)
% %clf
%  surf(grid.x,grid.y,flow.h)
%pcolor(x,y,Phi(:,:,end)')
%hold on
% %quiver3(X,Y,Phi(:,:,end)',qx,qy,zeros(size(qx)))
% quiver(X,Y,qx,qy)
% hold off
% xlabel('x [m]')
% ylabel('y [m]')
% colorbar
% title({'Girinskij potential \Phi',strcat('\rm{t = }',num2str(t(n)),' [s]')})

end

 figure(1)
% %clf
 surf(grid.x,grid.y,flow.zb,grid.x,grid.y,flow.h)
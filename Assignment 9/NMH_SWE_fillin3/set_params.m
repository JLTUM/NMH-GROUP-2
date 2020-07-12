function [grid, run, constants] = set_params(infilename)

    % Gravity [m/s^2]
    constants.g = 9.81;
    
    % number of ghost cell
    grid.NGHOST = 1;

    % reading parameters from input file and parsing to strucs

    load(infilename);
    
    grid.nx = nx;         
    grid.ny = ny;
    grid.xmax = xmax;       
    grid.xmin = xmin;
    grid.ymax = ymax;       
    grid.ymin = ymin;         
    
    run.dt = dt;       
    run.ntst = ntst;
     
    
    grid.nx = 250;
    grid.ny = 10;
    run.t = 6;
    run.dt = 0.005;
    run.ntst = run.t / run.dt;
    
    global infilename

% Case 1
if strcmp(infilename,'infile_2D_swe_damBreak.mat')
    
    grid.nx = 250;
    grid.ny = 10;
    run.t = 6;
    run.dt = 0.005;
    run.ntst = run.t / run.dt;
    
% Case 2
elseif strcmp(infilename,'infile_2D_swe_damBreak_V1.mat')
    
    grid.nx = 250;
    grid.ny = 10;
    run.t = 6;
    run.dt = 0.005;
    run.ntst = run.t / run.dt;

% Case 3
elseif strcmp(infilename,'infile_2D_swe_damBreak_V2.mat')
    
    grid.nx = 250;
    grid.ny = 10;
    run.t = 6;
    run.dt = 0.005;
    run.ntst = run.t / run.dt;
    
% Case 4
elseif strcmp(infilename,'infile_2D_swe_damBreak_V3.mat')
    
    grid.nx = 250;
    grid.ny = 10;
    run.t = 6;
    run.dt = 0.005;
    run.ntst = run.t / run.dt;
    
   % Case 5
elseif strcmp(infilename,'infile_2D_swe_damBreak_V4.mat')
    
    grid.nx = 250;
    grid.ny = 10;
    run.t = 6;
    run.dt = 0.005;
    run.ntst = run.t / run.dt;
    
    % Case 6
elseif strcmp(infilename,'infile_2D_swe_damBreak_V5.mat')
    
    grid.nx = 250;
    grid.ny = 10;
    run.t = 6;
    run.dt = 0.005;
    run.ntst = run.t / run.dt;
    
elseif strcmp(infilename,'infile_2D_swe_damBreak_V6.mat')
    
    grid.nx = 400;
    grid.ny = 10;
    run.t = 6;
    run.dt = 0.001;
    run.ntst = run.t / run.dt;
end
end




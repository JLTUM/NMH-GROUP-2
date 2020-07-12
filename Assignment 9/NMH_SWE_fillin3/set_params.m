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
end




function [grid, flow] = set_params(grid, flow, infilename)

    % reading parameters from input file and parsing to strucs

    load(infilename);
    
    grid.dt = dt;       
    grid.ntst = ntst;
    grid.nx = nx;         
    grid.ny = ny;
    grid.nu = nu;
    grid.xmax = xmax;       
    grid.xmin = xmin;
    grid.ymax = ymax;       
    grid.ymin = ymin;         
    
    %**********************************************************************
    % to manually override automatic input, uncomment below
    % and reset params:
    %
    % grid.dt = ...
    grid.ntst = 50;
    %
    %**********************************************************************
    
end





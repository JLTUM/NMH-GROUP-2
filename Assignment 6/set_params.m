
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
    grid.dt = 1e-4;
%     grid.ntst = 10000;
    grid.nx = 50;         
    grid.ny = 50;
    grid.nu = 0;
    %
    %**********************************************************************
    
end



function [grid, flow] = set_params(grid, flow, infilename)

    % reading parameters from input file and parsing to strucs

    load(infilename);
    
    grid.dt = 7e-4;       
    grid.ntst = 400;
    grid.nx = 30;         
    grid.ny = 30;
    grid.nu = 1;
    grid.xmax = xmax;       
    grid.xmin = xmin;
    grid.ymax = ymax;       
    grid.ymin = ymin;         
    
    %**********************************************************************
    % to manually override automatic input, uncomment below
    % and reset params:
    %
    %grid.dt = 2e-2;
    %grid.ntst = 100;
    %grid.nu = 10;
    %
    %**********************************************************************
    
end






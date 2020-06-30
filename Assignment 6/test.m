clear
close all

grid.xmin = 0;
grid.xmax = 1;
grid.ymin = 0;
grid.ymax = 1;
grid.nx = 25;
grid.ny = 25;


% Compute grid spacing
grid.dx = (grid.xmax - grid.xmin) / grid.nx;
grid.dy = (grid.ymax - grid.ymin) / grid.ny;

% Create coordinate arrays (xmin-dx/2 and xmax+dx/2) for ghost cells
grid.x = (grid.xmin - grid.dx/2) : grid.dx : (grid.xmax + grid.dx/2);
grid.y = (grid.ymin - grid.dy/2) : grid.dy : (grid.ymax + grid.dy/2);
   


    ii = 0 : grid.nx+2;
    jj = 0 : grid.ny+2;

    f = 7;
    h0 = 1;

    xsq = ((ii - grid.nx/2)*f/grid.nx).*((ii - grid.nx/2)*f/grid.nx); 
    ysq = ((jj - grid.ny/2)*f/grid.ny).*((jj - grid.ny/2)*f/grid.ny);

    for i = 1 : grid.nx+2

        for j = 1 : grid.ny+2
            flow.h(i,j) = h0 + 0.1*exp(-xsq(i))*exp(-ysq(j));
        end
    end
        
    size(grid.x)
    size(grid.y)
    size(flow.h)
    
surf(grid.x,grid.y,flow.h)
    
    
    
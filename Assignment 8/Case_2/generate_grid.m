function [ grid ] = generate_grid( grid )

%GENERATE_GRID Creates a struct containing all grid information.

% Compute grid spacing
grid.dx = ( grid.xmax - grid.xmin ) / grid.nx;
grid.dy = ( grid.ymax - grid.ymin ) / grid.ny;

% Create coordinate arrays
grid.x = ( grid.xmin - ( grid.NGHOST - 0.5 ) * grid.dx ):grid.dx:( grid.xmax + ( grid.NGHOST - 0.5 ) * grid.dx );
grid.y = ( grid.ymin - ( grid.NGHOST - 0.5 ) * grid.dy ):grid.dy:( grid.ymax + ( grid.NGHOST - 0.5 ) * grid.dy );

end


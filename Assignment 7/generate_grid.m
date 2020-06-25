function [ grid ] = generate_grid( grid )

%GENERATE_GRID Creates a struct containing all grid information.

% Compute grid spacing
grid.dx = (grid.xmax-grid.xmin)/grid.nx;
grid.dy = (grid.ymax-grid.ymin)/grid.ny;

% Create coordinate arrays
grid.x = (grid.xmin-grid.dx/2):grid.dx:(grid.xmax+grid.dx/2);
grid.y = (grid.ymin-grid.dy/2):grid.dy:(grid.ymax+grid.dy/2);

end


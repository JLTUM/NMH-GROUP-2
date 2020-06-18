function [ grid ] = generate_grid(grid)
%GENERATE_GRID Creates a struct containing all grid information.


% Compute grid spacing
grid.dx = ( grid.xmax - grid.xmin ) / grid.nx;
grid.dy = ( grid.ymax - grid.ymin ) / grid.ny;

% Create coordinate arrays
grid.x = grid.xmin :grid.dx: grid.xmax;
grid.y = grid.ymin :grid.dy: grid.ymax;

% general initialisation: ip, im, jp, jm, atbounds 

for i = 1 : grid.nx

    grid.ip(i) = i+1;
    grid.im(i) = i-1;

end

for j = 1 : grid.ny

    grid.jp(j) = j+1;
    grid.jm(j) = j-1;
    grid.atbounds(j) = 0;

end

% boundaries of domain (in the genral case)
% specific definitions follow in the switch section

grid.ip(grid.nx) = 1;
grid.jp(grid.ny) = 1;
grid.im(1) = grid.nx;
grid.jm(1) = grid.ny;
    

end


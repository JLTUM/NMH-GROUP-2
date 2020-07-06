function [ flow ] = set_initial_condition( grid, flow )
%SET_INITIAL_CONDITION Set initial fields

% ---- Create fields for shallow water equations ----------------------
% Flow depth (cell-centred)
flow.h = ones( length(grid.x), length(grid.y) );

% Specific discharge in x-direction (cell-centred)
flow.hu = ones( length(grid.x), length(grid.y) );

% Specific discharge in y-direction (cell-centred)
flow.hv = zeros( length(grid.x), length(grid.y) );

% Strickler value (cell-centred)
flow.kst = zeros( length(grid.x) - 2 * grid.NGHOST, length(grid.y) - 2 * grid.NGHOST );

% Bottom elevation (cell-centred)
flow.zb = zeros( length(grid.x) , length(grid.y) );

% Strickler value (cell-centred)
flow.kst = ones( grid.nx, grid.ny );

% ---- Fields initialization ----------------------

flow.I_s = -0.07;

% Bottom elevation
flow.zb = repmat(flow.I_s .* grid.x',1,length(grid.y));


% Strickler
kst = 30;
flow.kst = kst * ones( grid.nx, grid.ny );
flow.kst = readmatrix('kst2.txt')';
flow.zb = readmatrix('zb2.txt')';
end

function [ flow ] = set_initial_condition( grid, flow )
%SET_INITIAL_CONDITION Set initial fields

% ---- Create fields for shallow water equations ----------------------
% Flow depth (cell-centred)
flow.h = zeros( length(grid.x), length(grid.y) );

% Specific discharge in x-direction (cell-centred)
flow.hu = zeros( length(grid.x), length(grid.y) );

% Specific discharge in y-direction (cell-centred)
flow.hv = zeros( length(grid.x), length(grid.y) );

% Strickler value (cell-centred)
flow.kst = zeros( length(grid.x) - 2 * grid.NGHOST, length(grid.y) - 2 * grid.NGHOST );

% Bottom elevation (cell-centred)
flow.zb = zeros( length(grid.x) , length(grid.y) );

% Strickler value (cell-centred)
flow.kst = ones( grid.nx, grid.ny );

% ---- Fields initialization ----------------------
%% Bottom elevation
%%% add zb via function or files 

% zb Slope
 %flow.zb = 0.5 * grid.x' * grid.y;
 %flow.zb = 0.0 * grid.x' * grid.y;
 
 % zb Sharp edge
 %flow.zb(:,1:1:12) = 0.4;

 % zb Constant
flow.zb = 0.0 * grid.x' * grid.y;

%% Water level 
% give a proper intial flow depth
% flow.h(:) = ???;
% Constant Initial Condition
 flow.h = ones(grid.nx+2,grid.ny+2);

%% Strickler
% Constant Strickler
kst = 30;
flow.kst = kst *ones( grid.nx, grid.ny );




end


function [ flow ] = set_initial_condition( grid, flow )
%SET_INITIAL_CONDITION Set initial fields

% ---- Create fields for shallow water equations ----------------------
% Flow depth (cell-centred)
flow.h = zeros( length(grid.x), length(grid.y) );

% Specific discharge in x-direction (cell-centred)
flow.hu = ones( length(grid.x), length(grid.y) );
% flow.hu = 
% Specific discharge in y-direction (cell-centred)
flow.hv = zeros( length(grid.x), length(grid.y) );

% Strickler value (cell-centred)
flow.kst = zeros( length(grid.x) - 2 * grid.NGHOST, length(grid.y) - 2 * grid.NGHOST );

% Bottom elevation (cell-centred)
flow.zb = zeros( length(grid.x) , length(grid.y) );

% Strickler value (cell-centred)
 flow.kst = ones( grid.nx, grid.ny );
 
 flow.I_S = -0.001;

% ---- Fields initialization ----------------------
%% Bottom elevation
%%% add zb via function or files 

% zb Slope
%  flow.zb = 0 * grid.x' * grid.y;

 flow.zb = repmat(-0.001 .* grid.x',1,length(grid.y));
%  flow.zb(1,:)=-flow.zb(1,:)
 

 % zb Sharp edge
 %flow.zb(:,1:1:12) = 0.4;

 % zb Constant
%flow.zb = 10 * grid.x' * grid.y;

%% Water level 
% give a proper intial flow depth
 %flow.h(:) = 1;
 %NWabfluss
  flow.h = repmat(-0.001 .* grid.x'+1.2,1,length(grid.y));
 % Water level is drawn from a lognormal distribution (must be positive)
h0 = 1;
dh0 = 0.2;

 %% random log initial condition
  %  flow.h = lognrnd( log(h0^2/sqrt(h0^2+dh0^2)), sqrt(log(1+dh0^2/h0^2)), grid.nx+2, grid.ny+2 );

 
% Constant Initial Condition
% flow.h = ones(grid.nx+2,grid.ny+2);

%% Strickler
% Constant Strickler
kst = 30;
flow.kst = kst *ones( grid.nx, grid.ny );




end


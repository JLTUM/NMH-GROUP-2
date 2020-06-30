
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

% ---- Fields initialization ----------------------
% Bottom elevation
% TODO TODO TODO TODO TODO TODO TODO (try with different bottom elevation)
flow.zb = 0.05 * grid.x' * grid.y;

% keep this part -> "WALL" boundary condition
flow.zb(1,:) = flow.zb(2,:);
flow.zb(end,:) = flow.zb(end-1,:);
flow.zb(:,1) = flow.zb(:,2);
flow.zb(:,end) = flow.zb(:,end-1);

% Water level is drawn from a lognormal distribution (must be positive)
h0 = 1;
dh0 = 0.2;

% random log initial condition
     flow.h = lognrnd( log(h0^2/sqrt(h0^2+dh0^2)), sqrt(log(1+dh0^2/h0^2)), grid.nx+2, grid.ny+2 );

% % exp initial condition
%     ii = 0 : grid.nx+2;
%     jj = 0 : grid.ny+2;
%     f = 7;
%     xsq = ((ii - grid.nx/2) * f/grid.nx) .* ((ii - grid.nx/2) * f/grid.nx); 
%     ysq = ((jj - grid.ny/2) * f/grid.ny) .* ((jj - grid.ny/2) * f/grid.ny);
%     for i = 1 : grid.nx+2
%         for j = 1 : grid.ny+2
%             flow.h(i,j) = h0 + 1*exp(-xsq(i))*exp(-ysq(j));
%         end
%     end

% flow.h = ones(grid.nx+2,grid.ny+2);

% Constant Strickler
kst = 20;
flow.kst = kst *ones( grid.nx, grid.ny );

% ---- Create boundary conditions -----------------------------------------
bconds.bwest = {'WALL'};
bconds.beast = {'WALL'};
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

end


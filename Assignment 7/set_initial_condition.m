function [flow] = set_initial_condition( grid, flow)
%SET_INITIAL_CONDITION Set initial fields

    ii = 0 : grid.nx-1;
    jj = 0 : grid.ny/2-1;

    xsq = ((ii - grid.nx/2)*30/grid.nx).*((ii - grid.nx/2)*30/grid.nx); 
    ysq = ((jj - grid.ny/4)*30/grid.ny).*((jj - grid.ny/4)*30/grid.ny);

    for i = 1 : grid.nx

        for j = 1 : floor(grid.ny/2)
            flow.u(i,j) = exp(-xsq(i))*exp(-ysq(j));
        end

        for j = ceil(grid.ny/2) : grid.ny
            flow.u(i,j) = 0;
        end
    end

    for i = 1 : grid.nx
        for j = 1 : grid.ny

            ii = i-1;
            jj = j-1;

            flow.v(i,j) = 0;
            flow.tr_u(i,j) = -(jj-grid.ny/2)*(grid.dy * grid.ny)/grid.ny*100;
            flow.tr_v(i,j) =  (ii-grid.nx/2)*(grid.dx * grid.nx)/grid.nx*100;
        end
    end

end


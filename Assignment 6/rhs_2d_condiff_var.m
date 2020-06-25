function [rhsu, rhsv] = rhs_2d_condiff_var( grid, flow )
%CONVECTION2D_EQUATIONS Evaluate the time derivatives in the 2D convection
%equations

%% References
% We use a Finite-Difference scheme 

%% Preallocate for speed

% allocation for speed up

rhsu = zeros(grid.nx, grid.ny);
rhsv = zeros(grid.nx, grid.ny);
    
    % computation of right hand side for both equations

    for i = 1 : grid.nx
        for j = 1 : grid.ny
            
            ad = advection_u_lin2_var(grid, flow, i, j);
            diff = diffusion_u(grid, flow, i, j);
            rhsu(i,j) = -ad + diff;
            
            ad = advection_v_lin2_var(grid, flow, i, j);
            diff = diffusion_v(grid, flow, i, j);
            rhsv(i,j) = -ad + diff;
            
        end
    end
end



% additional functions
%**************************************************************************

function [adv_u] = advection_u_lin2_var(grid, flow, i, j)

    % for convenience

    ip = grid.ip(i);
    im = grid.im(i);
    jp = grid.jp(j);
    jm = grid.jm(j);
    dx = grid.dx;
    dy = grid.dy;
    
    % advection term 
    
    adv_u = (flow.tr_u(i,j) * (flow.u(ip, j) - flow.u(im, j))/(2*dx)) + ...
            (flow.tr_v(i,j) * (flow.u(i, jp) - flow.u(i, jm))/(2*dy));
    

end

function [adv_v] = advection_v_lin2_var(grid, flow, i, j)

    % for convenience

    ip = grid.ip(i);
    im = grid.im(i);
    jp = grid.jp(j);
    jm = grid.jm(j);
    dx = grid.dx;
    dy = grid.dy;
    
    % advection term 
    
    adv_v = (flow.tr_u(i,j) * (flow.v(ip, j) - flow.v(im, j))/(2*dx)) + ...
            (flow.tr_v(i,j) * (flow.v(i, jp) - flow.v(i, jm))/(2*dy));

end

function [diff_u] = diffusion_u(grid, flow, i, j)

    % returns the diffusive term for u at {x(i),y(j)}, i.e. the Laplacian
    % multiplied with the viscosity
    
    ip = grid.ip(i);
    im = grid.im(i);
    jp = grid.jp(j);
    jm = grid.jm(j);
    dx = grid.dx;
    dy = grid.dy;
    nu = grid.nu;
    
    diff_u = nu*((flow.u(im,j) - 2*flow.u(i,j) + flow.u(ip,j))/(dx^2) + ...
                 (flow.u(i,jm) - 2*flow.u(i,j) + flow.u(i,jp))/(dy^2));      
             
end

function [diff_v] = diffusion_v(grid, flow, i, j)

    % returns the diffusive term for v at {x(i),y(j)}, i.e. the Laplacian
    % multiplied with the viscosity

    ip = grid.ip(i);
    im = grid.im(i);
    jp = grid.jp(j);
    jm = grid.jm(j);
    dx = grid.dx;
    dy = grid.dy;
    nu = grid.nu;
    
    diff_v = nu*((flow.v(im,j) - 2*flow.v(i,j) + flow.v(ip,j))/dx^2 + ...
                 (flow.v(i,jm) - 2*flow.v(i,j) + flow.v(i,jp))/dy^2);
             
end





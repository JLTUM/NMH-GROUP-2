function [flow] = time_step_rk(grid, flow,type)
    %type 0 is runge kutta and type 1 is euler
    % timintegration follows here
    % be aware of the two methods, use commenting to choose
    % loop over timesteps 
    for itst = 1 : grid.ntst
        [flow.rhsu, flow.rhsv] = rhs_2d_condiff_var(grid, flow);
        disp(['time step ', num2str(itst)]);
        if type == 1
        
        [flow] = euler_2d_condiff_var(grid, flow);
        
        end
        
        if type == 0
            
        [flow] = rk_2d_condiff_var(grid, flow);
        
        end
        
        Plot_data
        
    end
end



function [flow] = rk_2d_condiff_var(grid, flow)

    % Implement low-storage 3rd-order Runge-Kutta time integration method
    % by filling the indicated missing part

    % 1st Runge-Kutta step
	
    % note: original RHS is computed outside of this function, and stored
    % in flow.rhsu, rhsv	
    
    for i = 1 : grid.nx
        for j = 1 : grid.ny
            
            S1_u(i,j) = grid.dt * flow.rhsu(i,j);
            S1_v(i,j) = grid.dt * flow.rhsv(i,j);
            flow.u(i,j) = flow.u(i,j) + 1/3 * S1_u(i,j);
            flow.v(i,j) = flow.v(i,j) + 1/3 * S1_v(i,j);
            
        end
    end
    
    % 2nd Runge-Kutta step (note: "secrhs" stands for second RHS)
    
    [flow.secrhsu, flow.secrhsv] = rhs_2d_condiff_var(grid, flow);
    
    for i = 1 : grid.nx
        for j = 1 : grid.ny
            
            %flow.secrhsu(i,j) = flow.secrhsu(i,j) + ???
            %flow.secrhsv(i,j) = flow.secrhsv(i,j) + ???
            
            S2_u(i,j) = grid.dt * flow.secrhsu(i,j) - 5/9 * S1_u(i,j);
            S2_v(i,j) = grid.dt * flow.secrhsv(i,j) - 5/9 * S1_v(i,j);
            flow.u(i,j) = flow.u(i,j) + 15/16 * S2_u(i,j);
            flow.v(i,j) = flow.v(i,j) + 15/16 * S2_v(i,j);
        
        end
    end
    
    % 3rd Runge-Kutta step
    
    [flow.rhsu, flow.rhsv] = rhs_2d_condiff_var(grid, flow);
    
    for i = 1 : grid.nx
        for j = 1 : grid.ny
            
            %flow.rhsu(i,j) = flow.rhsu(i,j) + ???
            %flow.rhsv(i,j) = flow.rhsv(i,j) + ???
            
            S3_u(i,j) = grid.dt * flow.rhsu(i,j) - 153/128 * S2_u(i,j);
            S3_v(i,j) = grid.dt * flow.rhsv(i,j) - 153/128 * S2_v(i,j);
            flow.u(i,j) = flow.u(i,j) + 8/15 * S3_u(i,j);
            flow.v(i,j) = flow.v(i,j) + 8/15 * S3_v(i,j);
            
        end
    end

end


function [flow] = euler_2d_condiff_var(grid, flow)

    for i = 1 : grid.nx
        for j = 1 : grid.ny
            
            flow.u(i,j) = flow.u(i,j) + grid.dt * flow.rhsu(i,j) * ...
                (1-grid.atbounds(j));
            
            flow.v(i,j) = flow.v(i,j) + grid.dt * flow.rhsv(i,j) * ...
                (1-grid.atbounds(j));
            
        end
    end
       
end


function [phi_out, phi_a_out, CFL]=conv_ee(U0,points,dt)

% Variable allocation
xend   = 2.0 * pi;
dx = xend / ( points - 1 );
x = 0.0 : dx : xend;
tsteps = 1000;
CFL = (U0*dt)/dx;
phi_out(1,:) = x; % x values for plotting
phi_a_out(1,:) = x; % x values for plotting

% Initial condition
 phi = sin(x);

% % Check initial field:
% plot(x, phi, 'r');
% title("initial field")
    
%% Explicit Euler

    for i = 1 : tsteps

        % Periodic boundary conditions at x=0:
        phinew(1) = phi(points-1);

        % Loop over grid points in space:
        for j = 2 : points - 1

            phinew(j) = phi(j) - U0*((phi(j+1)-phi(j-1))/(2*dx))*dt;

        end

        % Periodic boundary conditions at x=2*pi:
        phinew(points) = phi(2);

        % Write new field back to old field:
        phi = phinew;

        % Analytical solution
        phi_a = sin(x-U0*(i*dt));

        % write solution to matrix for plot every second time step
        if ~mod(i,2) == 1
            phi_out(end+1,:) = phi;
            phi_a_out(end+1,:) = phi_a;
            
        end
        
    end
    
end 



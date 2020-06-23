function [phi_out, phi_a_out, CFL]=conv_ie(U0,points,dt)

% Variable allocation
xend   = 2.0 * pi;
dx = xend / ( points - 1 );
x = 0.0 : dx : xend;
tsteps = 1000;
CFL = (U0*dt)/dx;
phi_out(1,:) = x; % x values for plotting
phi_a_out(1,:) = x; % x values for plotting
A = zeros(points,points);

% Initial Condition
phi = sin(x)';

% Check initial field:
% plot(x, phi, 'r');
% hold on;

%% Implicit Euler

% Compute coefficients of matrix A
a_w = -(U0*dt)/(2*dx);
a_p = 1;
a_e = (U0*dt)/(2*dx);

    for i = 1 : tsteps

        % Periodic boundary conditions at x=0:
        A(1,1) = a_p;
        A(1,2) = a_e;
        A(1,points-1) = a_w;
        b = phi;

        % Loop over grid points in space:
        for j = 2 : points - 1

            A(j,j) = a_p;
            A(j,j+1) = a_e;
            A(j,j-1) = a_w;

        end

        % Periodic boundary conditions at x=2*pi:
        A(points,points) = a_p;
        A(points,2) = a_e;
        A(points,points-1) = a_w;

        % Solve the linear system of equations
        phi = (A\b);

        % Analytical solution
        phi_a = sin(x - U0*(i*dt));

        % write solution to matrix for plot
        if ~mod(i,2) == 1

            phi_out(end+1,:) = phi;
            phi_a_out(end+1,:) = phi_a;

        end

    end

end


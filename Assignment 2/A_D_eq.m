%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 1
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

% ----------------------------------------------------------------------- %

% solution of

%       d phi          d^2 phi
%  -U0 ------- + Gamma ------- = 0
%        dx             dx^2

% u=0
% 0 <= x <= 2pi

function A_D_eq(U0,Gamma,points,scheme)


% Clear all variables and plots.
format long;
hold off;

% Set convection velocity
% U0 = [10];

% Set diffusivity
% Gamma = 1.0;

% Discrete spacing in space!
xend   = 2.0*pi;

% points = 10; 
dx     = xend/(points-1);
% Grid with x locations:
x = 0.0 : dx : xend;

% Initialization of field
phi = zeros(points,1);

% Initialization of matrix A
A = zeros(points,points);

% Initialization of vector b
b = zeros(points,1);

% Boundary condition
phi_0   = 0.0;
phi_end = 1.0;

% Loop over grid points in space
% note that boundary points are excluded
for i = 2 : points-1

    if scheme == "Central"
        % Central Sheme
        a_w = ((U0/2*dx)+(Gamma/dx^2));%dieses minus ist problematisch
        a_p = (-2*Gamma)/dx^2;
        a_e = -((U0/2*dx)-(Gamma/dx^2));
        
    elseif scheme == "Upwind"

        % Upwind
        a_w = ((U0/dx)+(Gamma/dx^2));
        a_p = ((-U0/dx)-((2*Gamma)/dx^2));
        a_e = (Gamma/dx^2);
    end
    
     
    % assign values to matrix A
    
    A(i,i-1) = a_w;
    A(i,i) = a_p;
    A(i,i+1) = a_e;
    
end

% Boundary conditions

% at i = 1
 A(1,1) = 1;
 b(1) = phi_0;

% at i = points
 A(points,points) = 1;
 b(points) = phi_end;

% Solution of the linear system
phi = A\b;

% Analytical solution

phi_analytic = (exp((U0.*x/Gamma))-1)/(exp((2*pi*U0)/Gamma)-1);

% error

% nn = points / 2;
%er = (phi_analytic(nn) - phi(nn)) / phi_analytic(nn);

% Plot the solution

plot(x,phi,'r', x,phi_analytic, 'g');

% Plot the error
% 
% figure
% loglog(points,er)
end
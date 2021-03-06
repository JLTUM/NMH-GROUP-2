%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 2
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

% ----------------------------------------------------------------------- %

% solution of

%       d phi          d^2 phi
%  -U0 ------- + Gamma ------- = 0
%        dx             dx^2

% u=0
% 0 <= x <= 2pi

function [er] = A_D_eq(U0,Gamma,points,nn,scheme)

global sp; % sp: subplot number 
global fig_A_D_eq;

% Clear all variables and plots.
format long;

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
        a_w = ((U0/(2*dx))+(Gamma/dx^2));
        a_p = -(2*Gamma)/dx^2;
        a_e = ((-U0/(2*dx))+(Gamma/dx^2));
        % assign values to matrix A
        A(i,i-1) = a_w;
        A(i,i) = a_p;
        A(i,i+1) = a_e;
    
    elseif scheme == "Upwind"
        
        % Upwind
        a_w = ((U0/dx)+(Gamma/dx^2));
        a_p = ((-U0/dx)-((2*Gamma)/dx^2));
        a_e = (Gamma/dx^2);
        % assign values to matrix A
        A(i,i-1) = a_w;
        A(i,i) = a_p;
        A(i,i+1) = a_e;
    
    elseif scheme == "Both"
        
        % Upwind 
        a_w_u = ((U0/dx)+(Gamma/dx^2));
        a_p_u = ((-U0/dx)-((2*Gamma)/dx^2));
        a_e_u = (Gamma/dx^2);
        % Central
        a_w_c = ((U0/(2*dx))+(Gamma/dx^2));
        a_p_c = -(2*Gamma)/dx^2;
        a_e_c = ((-U0/(2*dx))+(Gamma/dx^2));
        % assign values to matrix A
        A(i,i-1,1) = a_w_u;
        A(i,i,1) = a_p_u;
        A(i,i+1,1) = a_e_u;
        A(i,i-1,2) = a_w_c;
        A(i,i,2) = a_p_c;
        A(i,i+1,2) = a_e_c;
        
    end
end

% Boundary conditions

% at i = 1
 A(1,1,:) = 1;
 b(1) = phi_0;

% at i = points
 A(points,points,:) = 1;
 b(points) = phi_end;

% Analytical solution
phi_analytic = (exp((U0.*x/Gamma))-1)/(exp((2*pi*U0)/Gamma)-1);

% Solution of the linear system
if scheme == "Both"
    phi(:,1) = A(:,:,1)\b; % Upwind
    phi(:,2) = A(:,:,2)\b; % Central
else
    phi = A\b;
end

% Plot the solution

if sp == 1
    fig_A_D_eq = figure('units','normalized','outerposition',[0 0 1 1]);
end

set(0,'CurrentFigure',fig_A_D_eq)
subplot(2,2,sp)
hold on
if sp == 1
title('U: 10 Points: 5');
elseif sp == 2
title('U: -10 Points: 5');
elseif sp == 3
title('U: 10 Points: 51');
elseif sp == 4
title('U: -10 Points: 51');
end

if scheme == "Upwind"
    plot(x,phi,'r', x,phi_analytic, 'k');
    legend('Upwind','Analytic')
elseif scheme == "Central"
    plot(x,phi,'g', x,phi_analytic, 'k');
    legend('Central','Analytic')
elseif scheme == "Both"
    plot(x,phi(:,1),'r',x,phi(:,2),'g',x,phi_analytic,'k');
    legend('Upwind','Central','Analytic')
end

sp = sp +1;

end


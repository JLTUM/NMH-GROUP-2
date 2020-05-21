%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 3
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

% ----------------------------------------------------------------------- %

% This is a Matlab code to solve 1D steady advection-diffusion equation
% discritised by finite-volume schemes 
%
% differential form of advection-diffusion eqn.:
%
%       d phi          d^2 phi
%  -U0 ------- + Gamma ------- = 0
%        dx             dx^2
%
% You are asked to fill in the missing parts to complete the implementation.
% Missing parts are marked by ???
%
% Author: Yoshiyuki Sakai
% Email: yoshiyuki.sakai@tum.de
%%
function [phi_val]=A_D_eq_cells(U0,Gamma,cells,nn)

global sp; % sp: subplot number 

% Clear all variables and plots.
format long;

% Set up grid cells
xend = 2.0 * pi;
dx = xend/(cells); % dx for cells
dx_analytic = xend/(1000); % dx for analytic 

% Array of grid cell centre locations:
x = dx/2:dx:xend-(dx/2);
x_analytic = dx_analytic/2:dx_analytic:xend-(dx_analytic/2);

% Initialization of cell-averaged field
phi = zeros(cells,1);
% Initialization of matrix A
A = zeros(cells,cells);
% Initialization of vector b
b = zeros(cells,1);
% Boundary cell face values
phi_0   = 0.0;
phi_end = 1.0;

% Loop over grid cells
% NB: boundary cells are excluded
   
a_w = (U0/2)+Gamma/dx;
a_p = -(2*Gamma)/dx;
a_e = -(U0/2)+(Gamma/dx);
     
for i = 2 : cells-1

%     assign values to LHS matrix A
     A(i,i) = a_p;
     A(i,i-1) = a_w;
     A(i,i+1) = a_e;

%     assign values to RHS vector b
     b(i) = 0;
     
end

% Boundary conditions (Dirichlet at boundary cell faces)

% at i = 1
 A(1,1) = -(U0/2)-(3*Gamma)/dx; %Ap,BC left side
 A(1,2) = a_e;                  %Ae
 b(1) = -(U0+(2*Gamma)/dx)*phi_0; 

% at i = cells
 A(cells,cells) = (U0/2)-(3*Gamma)/dx; %Ap,BC right side
 A(cells,cells-1) = (U0/2)+(Gamma/dx); %Aw
 b(cells) = (U0-(2*Gamma)/dx)*phi_end;

% Solution of the linear system
phi = A\b;
phi = [phi_0; phi; phi_end]; % include boundary conditions 
x = [0 x xend]; % include boundary conditions 

% Compute analytical solution
phi_analytic = (exp((U0.*x_analytic/Gamma))-1)/(exp((2*pi*U0)/Gamma)-1); %%correct?
phi_analytic = [phi_0 phi_analytic phi_end]; % include boundary conditions 
x_analytic = [0 x_analytic xend]; % include boundary conditions 

%% Compute relative error Task 3
nn = ceil(cells / 2); 
%relative error
err_rel = (phi_analytic(nn) - phi(nn)) / phi_analytic(nn);

%% Compute mean error Task 4
err_mean =(sqrt((phi_analytic(nn) - phi(nn))^2))/phi_analytic(nn);


%% Plot the numerical and analytical solutions


subplot(3,2,sp)
hold on
plot(x, phi, '-.r', x_analytic, phi_analytic, '--k');
legend('Numerical','Analytic')
if sp == 1
title('U: 1 cells: 51');
elseif sp == 2
title('U: -1 cells: 51');
elseif sp == 3
title('U: 10 cells: 5');
elseif sp ==4
title('U: -10 cells: 5');
elseif sp == 5
title('U: 10 cells: 51');
elseif sp == 6
title('U: -10 cells: 51');
end

%% Plot relative and mean error



%set(0,'Currenigure',fig_A_D_eq)
% subplot(3,3,sp)
% hold on
% if sp == 1
% title('U: 10 Points: 5');
% elseif sp == 2
% title('U: -10 Points: 5');
% elseif sp == 3
% title('U: 10 Points: 51');
% elseif sp == 4
% title('U: -10 Points: 51');
% end
% 
% if scheme == "Upwind"
%     plot(x,phi,'r', x,phi_analytic, 'k');
%     legend('Upwind','Analytic')
% elseif scheme == "Central"
%     plot(x,phi,'g', x,phi_analytic, 'k');
%     legend('Central','Analytic')
% elseif scheme == "Both"
%     plot(x,phi(:,1),'r',x,phi(:,2),'g',x,phi_analytic,'k');
%     legend('Upwind','Central','Analytic')
% end
phi_val=[A,b]
sp = sp +1;


% Plot the error as function of dx in log-log scale
%???

end
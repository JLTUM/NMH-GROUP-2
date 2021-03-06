%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 3
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer
%% 

function [x, phi, x_analytic, phi_analytic, err_rel, err_mean] = A_D_FV(U0,Gamma,cells)

% Clear all variables and plots.
format long;

% Set up grid cells
xend = 2.0 * pi;
dx = xend/(cells); % dx for cells
dx_analytic = xend/(1001); % dx for analytic 

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
phi = [phi_0; phi; phi_end]'; % include boundary conditions 
x = [0 x xend]; % include boundary conditions 

% Compute analytical solution
x_analytic = [0 x_analytic xend]; % include boundary conditions 
phi_analytic = (exp((U0.*x_analytic/Gamma))-1)/(exp((2*pi*U0)/Gamma)-1);
phi_ex = (exp((U0.*x/Gamma))-1)/(exp((2*pi*U0)/Gamma)-1);

% Compute relative error
 nn = ceil(cells / 3 ) + 1; % value at pi
 err_rel = abs((phi_ex(nn) - phi(nn)) / phi_ex(nn));
 err_mean = sqrt(mean((phi_ex - phi).^2))/mean(phi_ex);
 
%  err_rel = abs((phi_ex - phi) ./ phi_ex);
%  plot(x,err_rel)
 
end
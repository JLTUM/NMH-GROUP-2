% Course: CFD Lab
% TU Muenchen, Summer term 2020
%Group 2
%Assignment 4 - Main code
%Andreas Mirlach, Julian Lenz, Faro Schäfer, Nick Pfeiffer
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi
%  ------- =  -U0 * -------
%    dt               dx
%
% 0 <= x <= 2pi
%
% periodic boundary condition
% phi(0) = phi(2pi)
%
% initial condition
% t = 0  ==>  phi = sin(x)
%
% Central difference scheme (CDS) for spatial discretization
% Explicit Euler scheme for time advancement
%
%--------------------------------------------------------------------------
%
clear ,clc, close all;
global sp;
sp = 1;
%% Task 1: explicit Euler time stepping

% Set convection velocity
U0 = 1.0;
% Discrete spacing in space
xend   = 2.0 * pi;
points = 40;
% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
dt     = 0.1;

[matrix]=conv_ee(U0,xend,points,tsteps,dt);

% %% homework 3 task 1
% 
% % U = 1
% U0 = 1;
% Gamma = 1;
% cells = [25 51 71 101 1001 10001]; 
% for i = 1:length(cells)
% [x, phi, x_analytic_1(1,:), phi_analytic_1(1,:), err_rel_1(i), err_mean_1(i)] = A_D_FV (U0,Gamma,cells(i));
% X_1{i,:} = x;
% PHI_1{i,:} = phi;
% CELLS_1(i) = cells(i);
% end
% 
% % U = 10
% U0 = 10;
% Gamma = 1;
% cells = [25 51 101]; 
% for i = 1:length(cells)
% [x, phi, x_analytic_2(1,:), phi_analytic_2(1,:), err_rel_2(i), err_mean_2(i)] = A_D_FV (U0,Gamma,cells(i));
% X_2{i,:} = x;
% PHI_2{i,:} = phi;
% CELLS_2(i) = cells(i);
% end
% 
% % U = -10
% U0 = -10;
% Gamma = 1;
% cells = [25 51 101]; 
% for i = 1:length(cells)
% [x, phi, x_analytic_3(1,:), phi_analytic_3(1,:), err_rel_3(i), err_mean_3(i)] = A_D_FV (U0,Gamma,cells(i));
% X_3{i,:} = x;
% PHI_3{i,:} = phi;
% CELLS_3(i) = cells(i);
% end
% 
% %% Plot
% 
% fig_AD_FV = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:3
%     
%     %U: +10
%     subplot(3,2,i+(i-1))
%     hold on
%     plot(X_2{i,:}, PHI_2{i,:}, '-.r', x_analytic_2(1,:), phi_analytic_2(1,:), '--k');
%     legend('Numerical','Analytic') 
%     xlabel('x')
%     xlabel('y')
%     if i == 1
%     title('U: 10 cells: 25');
%     elseif i == 2
%     title('U: 10 cells: 51');
%     elseif i == 3
%     title('U: 10 cells: 101');
%     end
%    
%     %U: -10
%     subplot(3,2,i+i)
%     plot(X_3{i,:}, PHI_3{i,:}, '-.r', x_analytic_3(1,:), phi_analytic_3(1,:), '--k');
%     legend('Numerical','Analytic') 
%     xlabel('x')
%     xlabel('y')
%     if i == 1
%     title('U: -10 cells: 25');
%     elseif i == 2
%     title('U: -10 cells: 51');
%     elseif i == 3
%     title('U: -10 cells: 101');
%     end
% 
% end
% 
% fig_error = figure('units','normalized','outerposition',[0 0 1 1]);
%     loglog((2*pi)./CELLS_1, err_rel_1,'-rx',(2*pi)./CELLS_1, err_mean_1,'-bo',...
%         (2*pi)./CELLS_1,((2*pi)./CELLS_1).^2,'k')
%     title('Error plot for U = 1')
%     legend('Relative error','RMS','Error order: 2')
%     xlabel('dx')
%     ylabel('RMS / Relative Error')
% 
% mkdir Plots_three
% print(fig_error,'-dpng',"Plots_three/Error.png",'-r150');
% print(fig_AD_FV,'-dpng',"Plots_three/Advection_diffusion_FV.png",'-r150');
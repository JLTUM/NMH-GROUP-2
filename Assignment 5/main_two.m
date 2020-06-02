%% Course: Numerical Methods
% TU Muenchen, Summer term 2020
% Group 2
% Assignment 5 - Main code
% Members: Andreas Mirlach, Julian Lenz, Faro Schäfer, Nick Pfeiffer
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi             d^2 phi
%  ------- =  -U0 * ------- + Gamma * ---------
%    dt               dx                dx^2 
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
% Crank Nicolson scheme for time advancement

clear, clc, close all;

% Preallocation of cells
CFL_cn = [];
Pe_cell_cn = [];
Phi_cn = cell(0);
Phi_a_cn = cell(0);

%% Input

% Case 1
U0 = 1;        points = 100;     dt = 0.1;    Gamma = 0;

disp('Calculating Case 1...')
[Phi_cn{end+1}, Phi_a_cn{end+1}, CFL_cn(end+1), Pe_cell_cn(end+1)] = conv_cn(Gamma,U0,points,dt);
disp([' CFL:  ',num2str(CFL_cn(end))])
disp([' Pe_cell:  ',num2str(Pe_cell_cn(end))])

% Case 2
U0 = 1;        points = 40;     dt = 0.1;    Gamma = 0;

disp('Calculating Case 2...')
[Phi_cn{end+1}, Phi_a_cn{end+1}, CFL_cn(end+1), Pe_cell_cn(end+1)] = conv_cn(Gamma,U0,points,dt);
disp([' CFL:  ',num2str(CFL_cn(end))])
disp([' Pe_cell:  ',num2str(Pe_cell_cn(end))])

% Case 3
U0 = 1;        points = 400;     dt = 1;    Gamma = 1;

disp('Calculating Case 3...')
[Phi_cn{end+1}, Phi_a_cn{end+1}, CFL_cn(end+1), Pe_cell_cn(end+1)] = conv_cn(Gamma,U0,points,dt);
disp([' CFL:  ',num2str(CFL_cn(end))])
disp([' Pe_cell:  ',num2str(Pe_cell_cn(end))])

% Case 4
U0 = 1;        points = 100;     dt = 0.1;    Gamma = 1;

disp('Calculating Case 4...')
[Phi_cn{end+1}, Phi_a_cn{end+1}, CFL_cn(end+1), Pe_cell_cn(end+1)] = conv_cn(Gamma,U0,points,dt);
disp([' CFL:  ',num2str(CFL_cn(end))])
disp([' Pe_cell:  ',num2str(Pe_cell_cn(end))])

%% Animation

%%% ------ uncomment everything in section Animation to see animation --- %%%

disp('Plotting...')

fig_animation = figure('units','normalized','outerposition',[0 0 1 1]); 

sp1 = subplot(221);
sp2 = subplot(222);
sp3 = subplot(223);
sp4 = subplot(224);

for i = 2:size(Phi_cn{1},1)-1

    % Plot Case 1
    set(gcf,'CurrentAxes',sp1)
        plot(Phi_cn{1}(1,:), Phi_cn{1}(i,:), Phi_a_cn{1}(1,:), Phi_a_cn{1}(i,:))

    % Plot Case 1
    set(gcf,'CurrentAxes',sp2)
        plot(Phi_cn{2}(1,:), Phi_cn{2}(i,:), Phi_a_cn{2}(1,:), Phi_a_cn{2}(i,:))
        
    % Plot Case 1
    set(gcf,'CurrentAxes',sp3)
        plot(Phi_cn{3}(1,:), Phi_cn{3}(i,:), Phi_a_cn{3}(1,:), Phi_a_cn{3}(i,:))
        
    % Plot Case 1
    set(gcf,'CurrentAxes',sp4)
        plot(Phi_cn{4}(1,:), Phi_cn{4}(i,:), Phi_a_cn{4}(1,:), Phi_a_cn{4}(i,:))
       
    drawnow limitrate
    pause(0.003)

end

%% Plot

fig_FP = figure('units','normalized','outerposition',[0 0 1 1]); 
    plot(Phi_cn{1}(1,:), Phi_cn{1}(400,:), Phi_a_cn{1}(1,:), Phi_a_cn{1}(400,:))
    legend('Crank Nicolson','Analytical')
    xlabel('x')
    ylabel('y')
    title('Floating Point Error')

fig_DispError = figure('units','normalized','outerposition',[0 0 1 1]); 
    plot(Phi_cn{2}(1,:), Phi_cn{2}(500,:), Phi_a_cn{2}(1,:), Phi_a_cn{2}(500,:))
    legend('Crank Nicolson','Analytical')
    xlabel('x')
    ylabel('y')
    title('Dispersive Error')
    
fig_Wiggles = figure('units','normalized','outerposition',[0 0 1 1]); 
    plot(Phi_cn{3}(1,:), Phi_cn{3}(400,:), Phi_a_cn{3}(1,:), Phi_a_cn{3}(400,:))
    legend('Crank Nicolson','Analytical')
    xlabel('x')
    ylabel('y')
    title('Wiggles')
    
% %% Print results
% 
% mkdir Plots_five
% print(fig_FP,'-dpng',"Plots_four/Damping_Explosion.png",'-r150');
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
Phi_cn = cell(0);
Phi_a_cn = cell(0);

%% Input

% Case 1
U0 = 0;        points = 40;     dt = 0.1;    Gamma = 1;

disp('Calculating Case 1...')
[Phi_cn{end+1}, Phi_a_cn{end+1}, CFL_cn(end+1)] = conv_cn(Gamma,U0,points,dt);
disp([' CFL:  ',num2str(CFL_cn(end))])

% Case 2
U0 = 1;        points = 40;     dt = 0.1;    Gamma = 0;

disp('Calculating Case 2...')
[Phi_cn{end+1}, Phi_a_cn{end+1}, CFL_cn(end+1)] = conv_cn(Gamma,U0,points,dt);
disp([' CFL:  ',num2str(CFL_cn(end))])

% Case 3
U0 = 1;        points = 100;     dt = 0.1;    Gamma = 1;

disp('Calculating Case 3...')
[Phi_cn{end+1}, Phi_a_cn{end+1}, CFL_cn(end+1)] = conv_cn(Gamma,U0,points,dt);
disp([' CFL:  ',num2str(CFL_cn(end))])

% Case 4
U0 = 1;        points = 4;     dt = 0.1;    Gamma = 1;

disp('Calculating Case 4...')
[Phi_cn{end+1}, Phi_a_cn{end+1}, CFL_cn(end+1)] = conv_cn(Gamma,U0,points,dt);
disp([' CFL:  ',num2str(CFL_cn(end))])

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
% 
% fig_Damping_Explosion = figure('units','normalized','outerposition',[0 0 1 1]); 
%     subplot(211)
%     plot(Phi_ee{1}(1,:), Phi_ee{1}(20,:), Phi_a_ee{1}(1,:), Phi_a_ee{1}(20,:))
%     legend('Explicit Euler','Analytical')
%     xlabel('x')
%     ylabel('y')
%     subplot(212)
%     plot(Phi_ie{1}(1,:), Phi_ie{1}(20,:), Phi_a_ie{1}(1,:), Phi_a_ie{1}(20,:))
%     legend('Implicit Euler','Analytical')
%     xlabel('x')
%     ylabel('y')
%     
% fig_Dispersive_Error = figure('units','normalized','outerposition',[0 0 1 1]);    
%     subplot(211)
%     plot(Phi_ee{2}(1,:), Phi_ee{2}(500,:), Phi_a_ee{2}(1,:), Phi_a_ee{2}(500,:))
%     legend('Explicit Euler','Analytical')
%     xlabel('x')
%     ylabel('y')
%     subplot(212)
%     plot(Phi_ie{2}(1,:), Phi_ie{2}(500,:), Phi_a_ie{2}(1,:), Phi_a_ie{2}(500,:))
%     legend('Implicit Euler','Analytical')
%     xlabel('x')
%     ylabel('y')
% 
% %% Print results
% 
% mkdir Plots_four
% print(fig_Damping_Explosion,'-dpng',"Plots_four/Damping_Explosion.png",'-r150');
% print(fig_Dispersive_Error,'-dpng',"Plots_four/Dispersive_error.png",'-r150');
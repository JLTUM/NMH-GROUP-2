%% Course: CFD Lab
% TU Muenchen, Summer term 2020
% Group 2
% Assignment 4 - Main code
% Members: Andreas Mirlach, Julian Lenz, Faro Schäfer, Nick Pfeiffer
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
% Explicit and Implicit Euler scheme for time advancement

clear, clc, close all;

% Preallocation of cells
CFL_ee = [];
Phi_ee = cell(0);
Phi_a_ee = cell(0);
CFL_ie = [];
Phi_ie = cell(0);
Phi_a_ie = cell(0);


%% Task 1: explicit Euler time stepping

% Case 1
U0 = 0.1;         points = 40;      dt = 0.01;

disp('Calculating Case 1...')
[Phi_ee{end+1}, Phi_a_ee{end+1}, CFL_ee(end+1)] = conv_ee(U0,points,dt);
[Phi_ie{end+1}, Phi_a_ie{end+1}, CFL_ie(end+1)] = conv_ie(U0,points,dt);
disp([' CFL:  ',num2str(CFL_ie(end))])

% Case 2
U0 = 1;        points = 100;     dt = 0.01;

disp('Calculating Case 2...')
[Phi_ee{end+1}, Phi_a_ee{end+1}, CFL_ee(end+1)] = conv_ee(U0,points,dt);
[Phi_ie{end+1}, Phi_a_ie{end+1}, CFL_ie(end+1)] = conv_ie(U0,points,dt);
disp([' CFL:  ',num2str(CFL_ie(end))])

% Case 3
U0 = 1;        points = 1000;     dt = 0.01;

disp('Calculating Case 3...')
[Phi_ee{end+1}, Phi_a_ee{end+1}, CFL_ee(end+1)] = conv_ee(U0,points,dt);
[Phi_ie{end+1}, Phi_a_ie{end+1}, CFL_ie(end+1)] = conv_ie(U0,points,dt);
disp([' CFL:  ',num2str(CFL_ie(end))])
 
% % Case 4
% U0 = 1;        points = 1000;      dt = 0.01;
% 
% disp('Calculating Case 4...')
% [Phi_ee{end+1}, Phi_a_ee{end+1}, CFL_ee(end+1)] = conv_ee(U0,points,dt);
% [Phi_ie{end+1}, Phi_a_ie{end+1}, CFL_ie(end+1)] = conv_ie(U0,points,dt);
% disp([' CFL:  ',num2str(CFL_ie(end))])


%% Animation

disp('Plotting...')
fig_animation = figure('units','normalized','outerposition',[0 0 1 1]); 

sp1 = subplot(4,2,1);
    title('U0=1 | dt=0.1 | points=1000')
sp2 = subplot(4,2,2);
    title('U0=1 | dt=0.1 | points=1000')
sp3 = subplot(4,2,3);
    title('U0=1 | dt=0.1 | points=1000')
sp4 = subplot(4,2,4);
    title('U0=1 | dt=0.1 | points=1000')
sp5 = subplot(4,2,5);
    title('U0=1 | dt=0.1 | points=1000')
sp6 = subplot(4,2,6);
    title('U0=1 | dt=0.1 | points=1000')
% sp7 = subplot(4,2,7);
%     title('U0=1 | dt=0.1 | points=1000')
% sp8 = subplot(4,2,8);
%     title('U0=1 | dt=0.1 | points=1000')
text(0.35,1.2,'column 1');    
sgtitle('Explicit Euler                                        Implicit Euler')

for i = 2:size(Phi_ee{1},1)-1

    % Plot Case 1
    set(gcf,'CurrentAxes',sp1)
        plot(Phi_ee{1}(1,:), Phi_ee{1}(i,:), Phi_a_ee{1}(1,:), Phi_a_ee{1}(i,:))
    set(gcf,'CurrentAxes',sp2)
        plot(Phi_ie{1}(1,:), Phi_ie{1}(i,:), Phi_a_ie{1}(1,:), Phi_a_ie{1}(i,:))
    
    % Plot Case 2
    set(gcf,'CurrentAxes',sp3)
        plot(Phi_ee{2}(1,:), Phi_ee{2}(i,:), Phi_a_ee{2}(1,:), Phi_a_ee{2}(i,:))
    set(gcf,'CurrentAxes',sp4)
        plot(Phi_ie{2}(1,:), Phi_ie{2}(i,:), Phi_a_ie{2}(1,:), Phi_a_ie{2}(i,:))
        
    % Plot Case 3
    set(gcf,'CurrentAxes',sp5)
        plot(Phi_ee{3}(1,:), Phi_ee{3}(i,:), Phi_a_ee{3}(1,:), Phi_a_ee{3}(i,:))
    set(gcf,'CurrentAxes',sp6)
        plot(Phi_ie{3}(1,:), Phi_ie{3}(i,:), Phi_a_ie{3}(1,:), Phi_a_ie{3}(i,:))
        
%     % Plot Case 4
%     set(gcf,'CurrentAxes',sp7)
%         plot(Phi_ee{4}(1,:), Phi_ee{4}(i,:), Phi_a_ee{4}(1,:), Phi_a_ee{4}(i,:))
%     set(gcf,'CurrentAxes',sp8)
%         plot(Phi_ie{4}(1,:), Phi_ie{4}(i,:), Phi_a_ie{4}(1,:), Phi_a_ie{4}(i,:))
       
    drawnow limitrate
    pause(0.003)

end
% keine ahung wieso die  muttergefickten titles nich in den Subplots sind ich versuchs seit 2 stunden
% animations nur für uns damit wir besser vergleichen können

%% Plot

% fig_EE_IE = figure('units','normalized','outerposition',[0 0 1 1]); 

% hier kommen die plots hin die wir dann ausprinten

%% Print results

% mkdir Plots_four
% print(fig_xxx,'-dpng',"Plots_three/Advection_diffusion_FV.png",'-r150');
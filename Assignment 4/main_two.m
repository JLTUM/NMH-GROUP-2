% Course: CFD Lab
% TU Muenchen, Summer term 2020
%Group 2
%Assignment 4 - Main code
%Andreas Mirlach, Julian Lenz, Faro Schäfer, Nick Pfeiffer

%% Course: CFD Lab
% TU Muenchen, Summer term 2020
% Group 2
% Assignment 4 - Main code
% Members: Andreas Mirlach, Julian Lenz, Faro Schäfer, Nick Pfeiffer

clear ,clc, close all;

%% Task 1: explicit Euler time stepping

U0 = 1;
xend   = 2.0 * pi;
points = 1000;
tsteps = 1000;
dt     = 0.01;

[Phi_ee, Phi_a_ee, error_ee, CFL_ee] = conv_ee(U0,xend,points,tsteps,dt);
[Phi_ie, Phi_a_ie, error_ie, CFL_ie] = conv_ie(U0,xend,points,tsteps,dt);

%% Plot Explicit Euler (Forward)

figure 
plot (error_ee(:,1), error_ee(:,2))
xlabel("timestep")
ylabel("RMS")

figure
hold on
subplot(2,2,1)
    plot(Phi_ee(1,:), Phi_ee(2,:), Phi_ee(1,:), Phi_a_ee(2,:))
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
subplot(2,2,2)
    plot(Phi_ee(1,:), Phi_ee(3,:), Phi_ee(1,:), Phi_a_ee(3,:))
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
subplot(2,2,3)
    plot(Phi_ee(1,:), Phi_ee(4,:), Phi_ee(1,:), Phi_a_ee(4,:))
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
subplot(2,2,4)
    plot(Phi_ee(1,:), Phi_ee(5,:), Phi_ee(1,:), Phi_a_ee(5,:))
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
sgtitle(['U0: ',num2str(U0),', points:',num2str(points)])


%% Plot Implicit Euler (Backward)

figure 
plot (error_ee(:,1), error_ee(:,2))
xlabel("timestep")
ylabel("RMS")

figure
hold on
subplot(2,2,1)
    plot(Phi_ee(1,:), Phi_ee(2,:), Phi_ee(1,:), Phi_a_ee(2,:))
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
subplot(2,2,2)
    plot(Phi_ee(1,:), Phi_ee(3,:), Phi_ee(1,:), Phi_a_ee(3,:))
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
subplot(2,2,3)
    plot(Phi_ee(1,:), Phi_ee(4,:), Phi_ee(1,:), Phi_a_ee(4,:))
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
subplot(2,2,4)
    plot(Phi_ee(1,:), Phi_ee(5,:), Phi_ee(1,:), Phi_a_ee(5,:))
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
sgtitle(['U0: ',num2str(U0),', points:',num2str(points)])

% mkdir Plots_four
% print(fig_error,'-dpng',"Plots_three/Error.png",'-r150');
% print(fig_AD_FV,'-dpng',"Plots_three/Advection_diffusion_FV.png",'-r150');
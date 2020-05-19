clear ,clc, close all;
global sp;
sp = 1;


%% create matrix of data for plots

% U = 1
U0 = 1;
Gamma = 1;
cells = [25 51 101]; %[25 51 71 101 1001 1001]
for i = 1:length(cells)
[x, phi, x_analytic_1(1,:), phi_analytic_1(1,:), err_rel_1(i), err_mean_1(i)] = A_D_FV (U0,Gamma,cells(i));
X_1{i,:} = x;
PHI_1{i,:} = phi;
CELLS_1(i) = cells(i);
end

% U = -1
U0 = -1;
Gamma = 1;
cells = [25 51 101]; %[25 51 71 101 1001 1001]
for i = 1:length(cells)
[x, phi, x_analytic_2(1,:), phi_analytic_2(1,:), err_rel_2(i), err_mean_2(i)] = A_D_FV (U0,Gamma,cells(i));
X_2{i,:} = x;
PHI_2{i,:} = phi;
CELLS_2(i) = cells(i);
end

% U = 10
U0 = 10;
Gamma = 1;
cells = [5 51];  % [25 51 101]
for i = 1:length(cells)
[x, phi, x_analytic_3(1,:), phi_analytic_3(1,:), err_rel_3(i), err_mean_3(i)] = A_D_FV (U0,Gamma,cells(i));
X_3{i,:} = x;
PHI_3{i,:} = phi;
CELLS_3(i) = cells(i);
end

% U = -10
U0 = -10;
Gamma = 1;
cells = [5 51]; % before [25 51 101]
for i = 1:length(cells)
[x, phi, x_analytic_4(1,:), phi_analytic_4(1,:), err_rel_4(i), err_mean_4(i)] = A_D_FV (U0,Gamma,cells(i));
X_4{i,:} = x;
PHI_4{i,:} = phi;
CELLS_4(i) = cells(i);
end

%% Plot - task 2

fig_AD_FV_task2 = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:2
    
    %U: +10
    subplot(2,2,i+(i-1))
    hold on
    plot(X_3{i,:}, PHI_3{i,:}, '-.r', x_analytic_3(1,:), phi_analytic_3(1,:), '--k');
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
    if i == 1
    title('U: 10 cells: 5');
    elseif i == 2
    title('U: 10 cells: 51');
    end
   
    %U: -10
    subplot(2,2,i+i)
    plot(X_4{i,:}, PHI_4{i,:}, '-.r', x_analytic_4(1,:), phi_analytic_4(1,:), '--k');
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
    if i == 1
    title('U: -10 cells: 5');
    elseif i == 2
    title('U: -10 cells: 51');
    end

end

%% Plot - task 3

fig_AD_FV = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:3
    
    %U: +1
    subplot(3,2,i+(i-1))
    hold on
    plot(X_1{i,:}, PHI_1{i,:}, '-.r', x_analytic_1(1,:), phi_analytic_1(1,:), '--k');
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
    if i == 1
    title('U: 1 cells: 25');
    elseif i == 2
    title('U: 1 cells: 51');
    elseif i == 3
    title('U: 1 cells: 101');
    end
   
    %U: -1
    subplot(3,2,i+i)
    plot(X_2{i,:}, PHI_2{i,:}, '-.r', x_analytic_2(1,:), phi_analytic_2(1,:), '--k');
    legend('Numerical','Analytic') 
    xlabel('x')
    xlabel('y')
    if i == 1
    title('U: -1 cells: 25');
    elseif i == 2
    title('U: -1 cells: 51');
    elseif i == 3
    title('U: -1 cells: 101');
    end

end

fig_error = figure('units','normalized','outerposition',[0 0 1 1]);
    loglog((2*pi)./CELLS_1, err_rel_1,'-rx',(2*pi)./CELLS_1, err_mean_1,'-bo',...
        (2*pi)./CELLS_1,((2*pi)./CELLS_1).^2,'k')
    title('Error plot for U = 1')
    legend('Relative error','RMS','Error order: 2')
    xlabel('dx')
    ylabel('RMS / Relative Error')

mkdir Plots_three
print(fig_error,'-dpng',"Plots_three/Error_task3.png",'-r150');
print(fig_AD_FV,'-dpng',"Plots_three/Advection_diffusion_FV_task3.png",'-r150');
print(fig_AD_FV_task2,'-dpng',"Plots_three/Advection_diffusion_FV_task2.png",'-r150');
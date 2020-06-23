
clear,clc, close all;
global sp;
sp = 1;

%% task 1

U0 = 1;
Gamma = 1;
cells = 51; 
[phi_val]=A_D_eq_cells(U0,Gamma,cells);

U0 = -1;
Gamma = 1;
cells = 51; 
[phi_val]=A_D_eq_cells(U0,Gamma,cells);


%% task 2 

U0 = 10;
Gamma = 1;
cells = 5; 
[phi_val]=A_D_eq_cells(U0,Gamma,cells);

U0 = -10;
Gamma = 1;
cells = 5; 
[phi_val]=A_D_eq_cells(U0,Gamma,cells);

U0 = 10;
Gamma = 1;
cells = 51; 
[phi_val]=A_D_eq_cells(U0,Gamma,cells);

U0 = -10;
Gamma = 1;
cells = 51; 
[phi_val]=A_D_eq_cells(U0,Gamma,cells);


%% task 3

cells = [25,51,101];
U0=1;
Gamma=1;
for k = 1 : length(cells)
    cell = cells(k);
    nn = ceil(cell/2);
    U0 = 10;
    [er_10(k,:)] = A_D_eq_cells_er(U0,Gamma,cell,nn);
end
% 
% % U-10
% for k = 1 : length(points)
%     point = points(k);
%     nn = ceil(point/2);
%     U0 = -10;
%     scheme = "Both";
%     [er_neg10(k,:)] = A_D_eq_er(U0,Gamma,point,nn,scheme);
% end

Plot
fig_error = figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(cells, er_10(:,1), cells, er_10(:,2))
    hold on 
    loglog(cells, er_neg10(:,1), cells, er_neg10(:,2))
    xlabel("number of points")
    ylabel("rel error")
    legend('U10 Upwind','U10 Central','-U10 Upwind','-U10 Central')
% 
% mkdir Plots_two
% print(fig_error,'-dpng',"Plots_two/Error10.png",'-r150');
% print(fig_A_D_eq,'-dpng',"Plots_two/Advection_diffusion.png",'-r150');

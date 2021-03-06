close all
clear all
global sp;
global fig_A_D_eq;
sp = 1;

%% task 2 3 4 5
U0 = 10;
Gamma = 1;
points = 5;
nn = ceil(points/2);
scheme = "Both"; % "Central" or "Upwind or "Both"
A_D_eq(U0,Gamma,points,nn,scheme);

U0 = -10;
Gamma = 1;
points = 5;
nn = ceil(points/2);
scheme = "Both"; % "Central" or "Upwind or "Both"
A_D_eq(U0,Gamma,points,nn,scheme);

U0 = 10;
Gamma = 1;
points = 51;
nn = ceil(points/2);
scheme = "Both"; % "Central" or "Upwind or "Both"
A_D_eq(U0,Gamma,points,nn,scheme);

U0 = -10;
Gamma = 1;
points = 51;
nn = ceil(points/2);
scheme = "Both"; % "Central" or "Upwind or "Both"
A_D_eq(U0,Gamma,points,nn,scheme);


%% task 7 

points = [51, 71, 101, 151, 201, 501, 1001];

% U10
for k = 1 : length(points)
    point = points(k);
    nn = ceil(point/2);
    U0 = 10;
    scheme = "Both";
    [er_10(k,:)] = A_D_eq_er(U0,Gamma,point,nn,scheme);
end

% U-10
for k = 1 : length(points)
    point = points(k);
    nn = ceil(point/2);
    U0 = -10;
    scheme = "Both";
    [er_neg10(k,:)] = A_D_eq_er(U0,Gamma,point,nn,scheme);
end

% Plot
fig_error = figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(points, er_10(:,1), points, er_10(:,2))
    hold on 
    loglog(points, er_neg10(:,1), points, er_neg10(:,2))
    xlabel("number of points")
    ylabel("rel error")
    legend('U10 Upwind','U10 Central','-U10 Upwind','-U10 Central')

mkdir Plots_two
print(fig_error,'-dpng',"Plots_two/Error10.png",'-r150');
print(fig_A_D_eq,'-dpng',"Plots_two/Advection_diffusion.png",'-r150');
    
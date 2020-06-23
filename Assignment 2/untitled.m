% close all
global sp;
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

points = [5, 51, 71, 101, 151, 201];

for k = 1 : length(points)
    point = points(k);
    nn = ceil(point/2);
    U0 = 10;
    scheme = "Both";
    [phi_U10{:,:,k}, er_10(k,:)] = A_D_eq(U0,Gamma,point,nn,scheme);
end

for k = 1 : length(points)
    point = points(k);
    nn = ceil(point/2);
    U0 = -10;
    scheme = "Both";
    [phi_Uneg10{:,:,k}, er_Uneg10(k,:)] = A_D_eq_er(U0,Gamma,point,nn,scheme);
end



%% Plot task 2,3,4


subplot(2,2,sp)
if scheme == "Upwind"
    plot(x,phi,'r', x,phi_analytic, 'k');
%     hold on
elseif scheme == "Central"
    plot(x,phi,'g', x,phi_analytic, 'k');
%     hold on
elseif scheme == "Both"
    plot(x,phi(:,1),'r',x,phi(:,2),'g',x,phi_analytic,'k');
%     hold on
end



figure
plot(points, er_10(:,1), points, er_10(:,2))
hold on 
plot(points, er_neg10(:,1), points, er_neg10(:,2))
legend('U10 Upwind','U10 Central','-U10 Upwind','-U10 Central')


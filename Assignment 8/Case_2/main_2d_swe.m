%**************************************************************************
% NMH Lab Summer Semester 2020, Assignment 8
%
% This code solves the 2D shallow water equations
%
% author: H. Zeng & L. Unglehrt
% July, 2020
%**************************************************************************
clear;
close all

%% Initialize simulation
% read infile 
global infilename
infilename = 'infile_2D_swe_channelFlow2.mat'; %% 1,2,3,4,5
fprintf('infilename is: %s\n', infilename)

% build structures 
[grid, run, constants, flow, bconds] = build_structs;
fprintf('struct built\n')

% fill some fields of 'grid' and 'flow' with data from infile
[grid, run, constants] = set_params(infilename);
fprintf('parameters set\n')

% Generate an equidistant grid 
[grid] = generate_grid(grid);    

% Set initial conditions 
run.t = 0;
% ---- Create boundary conditions -----------------------------------------
[ bconds ] = set_boundary_conditions();

[ flow ] = set_initial_condition( grid, flow);

%% Preallocation of variables

b = grid.ymax;
Q = bconds.huwest*b;

I_WSP = [];
I_S = [];
I_E = [];
deltaE = [];
ymin = [];
ymincor=0.415; %%case 2 depth on bump

I = -flow.I_S;
k_st = flow.kst(1,1);
N_M = NWV_muster(Q,b,0,I,k_st);

%% Preallocation of Plots

    fig_Hy = figure;
        hold on
        title('H-y Diagram')
        y = 0:0.01:4;
        H = y + ( Q )^2 ./ (b * y ).^2 ./ 2 ./ 9.81;
        plot(H,y);
        H = y;
        y_c = (((Q/b)^2)/09.81)^(1/3);
        yline(y_c);
        plot(H,y);
        xlim([0 4])
        xlabel('H') 
        ylabel('y')
    fig_Channeld = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
    %caculate height of bump 
    
%% Time integration

for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );
    
    % Berechnung NWV
    u = flow.hu(:,2) ./ flow.h(:,2);
    A = flow.h(:,2) .* grid.ymax;
    U = 2 .* flow.h(:,2) + grid.ymax;
    R = A ./ U;
    v_st = k_st*sqrt(I) * R.^(2/3);
    
    % Berechnung Fr und H
    Fr = u ./ (flow.h(:,2).*9.81).^(0.5);
    H = u.^(2)/9.81/2 + flow.h(:,2) + flow.zb(:,2);
    
    % Berechnung Gefällelinien
    I_WSP(end+1) = abs((flow.zb(2,2) + flow.h(2,2)) - ((flow.zb(end,2) + flow.h(end,2))) )  / grid.xmax;
    I_S(end+1) = abs((flow.zb(2,2)-flow.zb(end))) / grid.xmax;
    I_E(end+1) = abs(H(2) - H(end)) / grid.xmax;
    
    %Berechnung Energieverlust
    deltaE(end+1)=abs(H(2) - H(end));
    ymin(end+1)=flow.h(51,2);
    
%% Plot Results
    
    % Plot Channel Diagnosis
    set(0, 'CurrentFigure', fig_Channeld)
    subplot(131)
    plot(grid.x(2:end), v_st(2:end), grid.x(2:end), u(2:end), grid.x(2:end), flow.hu(2:end,2));
    legend('v_{st}','u','hu')
    title('Channel Velocity')
    
    subplot(132)
    plot(grid.x(2:end), H(2:end), grid.x(2:end), flow.h(2:end,2)+flow.zb(2:end,2),...
        grid.x(2:end), flow.h(2:end,2), grid.x(2:end),flow.zb (2:end,2));
    title('Channel Waterdepth / Energy / Discharge')
    legend('H','h+zb','h','zb','Location','northeast')
    
    hold off
    subplot(133)
    plot(grid.x(2:end),Fr(2:end));
    title('Channel Froude')
    legend('Fr','Location','northwest')

  
    % Plot Hy-Diagramm
    set(0, 'CurrentFigure', fig_Hy)
    hold on
    scatter(H(51),flow.h(51,2),'b.')
    scatter(H(5),flow.h(5,2),'r.')
    scatter(H(100),flow.h(100,2),'g.')
    scatter(N_M(1),N_M(2),'yx')
    hold off
      
    pause(0.0001)
    
    %plot Energieverlust over time
    figure(9)
    plot(1:itstep,deltaE,'r',1:itstep,ymin,'b',1:itstep,ymincor,'g')
    legend('deltaE','y(min)','y(min) correct')
 
    


end
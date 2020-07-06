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

N_M = NWV_muster(Q,b,0,-flow.I_S,30);

I = -flow.I_S;
k_st = flow.kst(1,1);
N_M = NWV_muster(Q,b,0,I,k_st);

%% Preallocation of Plots

    fig_Hy = figure;
    fig_Channeld = figure;
        title('H-y Diagram')
        y = 0:0.01:5;
        hold on
        H = y + ( Q )^2 ./ (b * y ).^2 ./ 2 ./ 9.81;
        plot(H,y);
        H = y;
        y_c = (((Q/b)^2)/09.81)^(1/3);
        yline(y_c);
        plot(H,y);
        xlim([0 5])
        xlabel('H') 
        ylabel('y')

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
    
%% Plot Results
    
    % Plot Channel Diagnosis
    set(0, 'CurrentFigure', fig_Channeld)
    subplot(2,2,1)
    plot(grid.x(2:end),v_st(2:end),grid.x(2:end),u(2:end));
    title('Channel Velocity')
 
    subplot(2,2,2)
    hold on
    plot(1:itstep,I_WSP,1:itstep,I_S,1:itstep,I_E)
    hold off
    title('Channel I_WSP')
    
    subplot(2,2,3)
    plot(grid.x,flow.h(:,2)+flow.zb(:,2),grid.x,H,grid.x,flow.hu,grid.x,flow.zb,'b');
    title('Channel Waterdepth / Energy / Discharge')
    legend('h','H','q','zb','Location','southwest')
    
    hold off
    subplot(2,2,4)
    plot(grid.x,Fr);
    title('Channel Froude')
    legend('Fr','Location','northwest')

  
    % Plot Hy-Diagramm
    set(0, 'CurrentFigure', fig_Hy)
    hold on
    scatter(flow.h(90,2),H(90),'b')
    scatter(flow.h(20,2),H(20),'r')
    scatter(N_M(1),N_M(2),'y')
    hold off
      
    pause(0.0001)
    


end


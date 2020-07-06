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



% ---- Setup of for time integration --------------------------------------
% Frequency of diagnostic output
%itdiag = 100;

%% Preallocation of variables
b = grid.ymax;
Q = bconds.huwest*b;


%% Time integration

N_M = NWV_muster(Q,b,0,-flow.I_S,30);

%flow.h(:) = 1%N_M(1);
%flow.hu(:) = 0%N_M(3)*N_M(1);
I_WSP = [];
I_S = [];
I_E = [];
for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );
    
    I = -flow.I_S;
    k_st = flow.kst(1,1);
    N_M = NWV_muster(Q,b,0,I,k_st);
    v_h = flow.hu(:,2) ./ flow.h(:,2);
    A = flow.h(:,2) .* grid.ymax;
    U = 2 .* flow.h(:,2) + grid.ymax;
    R = A ./ U;
    v_st = k_st*sqrt(I)*R.^(2/3);
    
    Fr = v_h./(flow.h(:,2).*9.81).^(0.5);
    energy = v_h.^(2)/9.81/2+flow.h(:,2)+flow.zb(:,2);
    
    
    

    I_WSP(end+1) = abs((flow.zb(grid.NGHOST+1,grid.NGHOST+1) + flow.h(grid.NGHOST+1,grid.NGHOST+1)) -...
        ((flow.zb(end,grid.NGHOST+1) + flow.h(end,grid.NGHOST+1))) )  / grid.xmax;
    
    I_S(end+1) = abs((flow.zb(grid.NGHOST+1,grid.NGHOST+1)-flow.zb(end))) / grid.xmax;
   
    I_E(end+1) = abs(energy(grid.NGHOST+1) - energy(end)) / grid.xmax;
    

    %% Calculate NWV   
    

    
    figure(1)
    
    subplot(2,2,1)
    plot(grid.x(2:end),v_st(2:end),grid.x(2:end),v_h(2:end));
    title('Channel Velocity')
 
    subplot(2,2,2)
    plot(1:itstep,I_WSP,1:itstep,I_S,1:itstep,I_E)
    title('Channel I_WSP')
    
    subplot(2,2,3)
    plot(grid.x,flow.h(:,2)+flow.zb(:,2),grid.x,energy,grid.x,flow.hu,grid.x,flow.zb,'b');
    title('Channel Waterdepth / Energy / Discharge')
    legend('h','H','q','zb','Location','southwest')
    
    subplot(2,2,4)
    plot(grid.x,Fr);
    title('Channel Froude')
    legend('Fr','Location','northwest')
  
    figure(2)
    title('H-y Diagram')
    y = 0:0.01:5;
    %H = zeros(1,91);
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
    scatter(flow.h(90,2),energy(90))
    scatter(flow.h(20,2),energy(20))
    scatter(N_M(1),N_M(2),'b')
    hold off
    
    
 

%% Plot Results

%     v_n(end+1)=N_M(3);
%     % Channel Diagnosis
%     figure(1);
%     subplot(2,2,1);
%     plot(1:itstep,I_S,'-b', 1:itstep,I_WSP,'-r', 1:itstep, v_st,'-g')
%     legend('I_S','I_{WSP}','v_{st}','Location','northwest')
%     subplot(2,2,2);
%     plot(grid.x,flow.h(:,2)+flow.zb(:,2))
%     title('Channel diagnosis')
%     
%     subplot(2,2,3);
%     plot(grid.x,flow.hu(:,2)./ flow.h(:,2))
%     yline(N_M(3))
%     title('Channel Velocity Test')
%     
%     subplot(2,2,4);
%     energy = ( flow.hu(:,2)./flow.h(:,2) ).^2 / (2*9.81) + flow.h(:,2);
%     set(0, 'CurrentFigure', fig_hu)
%     plot(grid.x,energy)
%     yline(N_M(2))
%     title('Channel Energy Test')
%     
    pause(0.0001)
    


end


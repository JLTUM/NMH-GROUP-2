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
infilename = 'infile_2D_swe_channelFlow1.mat'; %% 1,2,3,4,5
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
[ flow ] = set_initial_condition( grid, flow );

% ---- Create boundary conditions -----------------------------------------
[ bconds ] = set_boundary_conditions();

% ---- Setup of for time integration --------------------------------------
% Frequency of diagnostic output
itdiag = 100;

%% Preallocation of variables

R_hyd = [];
v_st = [];
I_WSP = [];
I_S = [];
v_n = [];
% TODO: Froude number
% TODO: u_mean
b = grid.ymax;
Q = bconds.hwest*bconds.huwest*b;


%% Time integration
fprintf('start time integration\n')

for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );
    
    
    v_h = flow.hu(:,2) ./ flow.h(:,2);

    I = -flow.I_s;
    k_st = flow.kst(1,1);
    N_M = NWV_muster(Q,b,0,I,k_st);

    A = flow.h(:,2) .* grid.ymax;
    U = 2 .* flow.h(:,2) + grid.ymax;
    R = A ./ U;
    v_st = k_st*sqrt(abs(flow.I_s))*R.^(2/3);
    
    I_WSP = mean(mean(v_h))/ k_st / mean(mean(R))^(4/3);
    
    figure(1)
    hold on
    plot(grid.x,v_st,grid.x,v_h);
    
    yline(N_M(3));
    
    title('Channel Velocity')
    hold off
    
    
    
    %v_st(end+1) = min(min(flow.kst)) * sqrt(abs(flow.I_s)) *  R_hyd(end)^(2/3);
    %I_WSP(end+1) = abs((min(min(flow.h(2:end,:)))-max(max((flow.h(2:end,:))))))/grid.xmax;
    %I_S(end+1) = abs((min(min(flow.zb))-max(max(flow.zb))))/grid.xmax;

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
    pause(0.01)
    


end


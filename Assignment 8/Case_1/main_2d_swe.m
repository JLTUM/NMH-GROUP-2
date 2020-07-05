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
% ---- Create boundary conditions -----------------------------------------
[ bconds ] = set_boundary_conditions();

[ flow ] = set_initial_condition( grid, flow);



% ---- Setup of for time integration --------------------------------------
% Frequency of diagnostic output
%itdiag = 100;

%% Preallocation of variables
b = grid.ymax;
Q = bconds.huwest*b;
I = 0.001;
k_st = 30;

%% Time integration

N_M = NWV_muster(Q,b,0,I,k_st);

flow.h(:) = 10;
%flow.hu(:) = N_M(3);

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
    v_st = k_st*sqrt(I)*R.^(2/3);
    Fr = v_h./(flow.h(:,2).*9.81).^(0.5);
    %% Calculate NWV
    
    for i = 1:length(grid.x)-1
    I_WSP(i) = ( (v_h(i)+v_h(i+1)) /2 )/ k_st / ( (R(i)+R(i+1))/2 )^(4/3);
    end
    
    energy = v_h.^(2)/9.81/2+flow.h(:,2)+flow.zb(:,2);
    
    figure(1)
    
    subplot(2,2,1)
    plot(grid.x,v_st,grid.x,v_h);
    title('Channel Velocity')
 
    subplot(2,2,2)
    plot(grid.x,-I*grid.x,grid.x(2:end),I_WSP)
    title('Channel I_WSP')
    
    subplot(2,2,3)
    plot(grid.x,flow.h(:,2),grid.x,energy);
    title('Channel Waterdepth / Energy')
    legend('Depth','Energy','Location','northwest')
    
    subplot(2,2,4)
    plot(grid.x,Fr);
    title('Channel Waterdepth / Energy')
    legend('Fr','Location','northwest')
    
    
 

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


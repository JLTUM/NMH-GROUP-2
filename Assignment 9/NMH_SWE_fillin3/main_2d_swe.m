%**************************************************************************
% NMH Lab Summer Semester 2020, Assignment 9
%
% This code solves the 2D shallow water equations
%
% author: H. Zeng & L. Unglehrt
% July, 2020
%**************************************************************************
clear;
close all

%% Initialize simulation
global infilename
infilename = 'infile_2D_swe_damBreak.mat'; 
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

%% ---- Create boundary conditions -----------------------------------------
[ bconds ] = set_boundary_conditions();

% ---- Setup of for time integration --------------------------------------
% Frequency of diagnostic output
% itdiag = 


%% Time integration
itdiag=1;
fprintf('start time integration\n')
for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );

    % Diagonostic output
if mod(itstep, itdiag) == 0
        % CFL number
        fprintf('%d : CFL number:         %e\n', itstep, ...
            compute_CFL_number(constants, grid, run.dt, flow.h, flow.hu, flow.hv));
        
%% Plot results
%% Plot Results
% disp(itstep)
if ~mod(itstep,1)
    
   

    % Plot Channel Diagnosis
%    set(0, 'CurrentFigure', fig_Channeld)
%     subplot(131)
%     plot(grid.x(2:end), v_st(2:end), grid.x(2:end), u(2:end), grid.x(2:end), flow.hu(2:end,2));
%     legend('v_{st}','u','hu')
%     title('Channel Velocity')

%     subplot(131)
%     title('Channel Diagnosis')
%     yyaxis left
%     plot(1:itstep, Q_channeld)
%     ylabel('[m³/s]')
%     yyaxis right
%     plot(1:itstep, U_channeld)
%     ylabel('[m/s]')
%     xlabel('timestep')
%     legend('Q','U_b','Location','southwest')
%     hold off
    
%     subplot(1)
%     plot(grid.x(2:end), flow.h(2:end,2)+flow.zb(2:end,2),...
%          grid.x(2:end),flow.h (2:end,2), grid.x(2:end),flow.zb (2:end,2));
%     title('Head and Waterlevel')
%     legend('h+zb','h','zb','Location','southwest')
%     xlabel('x')
%     ylabel('[m]')
%     hold off
    
    surf(grid.x,grid.y,flow.h','FaceColor','b')
    
        xlabel('x','Fontsize',15)
        ylabel('y','Fontsize',15)
        zlabel('h','Fontsize',15)
%         a = get(gca,'XTickLabel');
%         set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')
%         set(gca,'YTickLabel',a,'fontsize',15,'FontWeight','bold')
%         set(gca,'ZTickLabel',a,'fontsize',15,'FontWeight','bold')
%         zlim([-1 2])
        title(['n=',num2str(itstep)])
        pause(0.05)
        %hold off  
        
%     subplot(133)
%     title('Froud and specific Discharge')
%     yyaxis left 
%     plot(grid.x(2:end), Fr(2:end));
%     ylabel('Froude [-]')
%     yyaxis right
%     plot(grid.x(2:end), flow.hu(2:end,2))
%     ylabel('[m²/s]')
%     xlabel('x')
%     legend('Fr','hu','Location','northwest')
%     hold off

    

    
end

    end
end

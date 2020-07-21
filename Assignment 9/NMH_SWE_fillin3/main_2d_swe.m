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
% read infile 
global infilename 
infilename = "infile_2D_swe_damBreak_V1.mat"; %%V1,V2;V3,V4,V5
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
[ flow1,flow2,flow3 ] = set_initial_condition( grid, flow );

% Create boundary conditions 
[ bconds ] = set_boundary_conditions();

%% Preallocation of variables

itdiag = 10;

U_channeld = [];
Q_channeld = [];

Frame = 1;
Video(run.ntst/itdiag) = struct('cdata',[],'colormap',[]);

VidObj = VideoWriter('Video');  % Name of Video
open(VidObj)

%% Preallocation of Plots

 fig_Channeld = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
 Axes = gcf;
 xlim([1 10]);
 ylim([-0.1 1.2]);
%  ax = gca;
%  ax.NextPlot = 'replaceChildren';

%% Time integration
dir= 'Plot_nine_'+infilename
mkdir (sprintf(dir)) 


     set(0, 'CurrentFigure', fig_Channeld)
        plot( grid.x(2:end),flow1.h (2:end,2),grid.x(2:end),flow2.h (2:end,2),grid.x(2:end),flow3.h (2:end,2));
        title(['initial setup'])
        %legend('h=1m','h=0.8m','h=0.4m','Location','north') %V1
        %legend('h(kst=10)','h(kst=30)','h(kst=110)','Location','north') %V2
        %legend('h0=0.05m','h0=0.1m','h0=0.2m','Location','north') %V1
        xlabel('x')
        ylabel('[m]')
        xlim([1 10]);
        ylim([-0.1 1]);
        hold off
        pause(0.001)
        
        
        
        print(fig_Channeld,'-dpng',fullfile(dir,"initial setup"),'-r150');
        

fprintf('start time integration\n')
for itstep = 1:run.ntst
    
    [ run, flow1 ] = time_step_rk1( itstep==1, constants, grid, run, ...
        flow1, bconds );

    [ run, flow2 ] = time_step_rk2( itstep==1, constants, grid, run, ...
        flow2, bconds );
    
    [ run, flow3 ] = time_step_rk3( itstep==1, constants, grid, run, ...
        flow3, bconds );
    % Calculation Channel Diagnosis
%     U_channeld(end+1) = mean(flow.hu(2:end-1,2) ./ flow.h(2:end-1,2));
%     Q_channeld(end+1) = mean(flow.hu(2:end-1,2) * grid.ymax);
    
   
    if mod(itstep, itdiag) == 0

    % Diagonostic output
        
        % Energy and flow head
        u = flow1.hu(:,2) ./ flow1.h(:,2);                            % u velocity
        v = flow1.hv(:,2) ./ flow1.h(:,2);                            % v velocity
        U = sqrt (u.^2 + v.^2);                                     % Velocity Vector
        H = U.^2 / (2*constants.g) + flow1.h(:,2) + flow1.zb(:,2);    % energy head
        v_h = U.^2 / (2*constants.g);                               % velocity head
        c = sqrt(constants.g .* flow1.h(:,2));
        Fr=u.*(c.^-1);
               
        % CFL number
        fprintf('%d : CFL number:         %e\n', itstep*run.dt, ...
            compute_CFL_number(constants, grid, run.dt, flow1.h, flow1.hu, flow1.hv));
        
     % Plot results
     H1=flow1.h;%/1;
     H2=flow2.h;%/0.8;
     H3=flow3.h;%/0.4;
    
     
% %    for parallel plotting  

%               set(0, 'CurrentFigure', fig_Channeld)
%         plot( grid.x(2:end),H1 (2:end,2),grid.x(2:end),H2 (2:end,2),grid.x(2:end),H3 (2:end,2));
%         title(['time = ',num2str(itstep * run.dt),'s'])
%         %legend('h=1m','h=0.8m','h=0.4m','Location','north') %V1
%         %legend('h(kst=10)','h(kst=30)','h(kst=110)','Location','north') %V2
%         legend('h0=0.05m','h0=0.1m','h0=0.2m','Location','north') %V1
%         xlabel('x')
%         ylabel('[m]')
%         xlim([1 10]);
%         ylim([-0.1 1]);
%         hold off
%         pause(0.001)
        
        set(0, 'CurrentFigure', fig_Channeld)
        plot( grid.x(2:end),flow1.h (2:end,2),grid.x(2:end),u(2:end),grid.x(2:end),c(2:end));
        title(['time = ',num2str(itstep * run.dt),'s'])
        %legend('h=1m','h=0.8m','h=0.4m','Location','north') %V1
        %legend('h(kst=10)','h(kst=30)','h(kst=110)','Location','north') %V2
        %legend('h0=0.05m','h0=0.1m','h0=0.2m','Location','north') %V1
        legend('h','u','c','Location','north') %V1
        xlabel('x')
        ylabel('[m]')
        xlim([1 10]);
        ylim([-0.5 5]);
        hold off
        pause(0.001)
     
     
%         set(0, 'CurrentFigure', fig_Channeld)
%         plot(grid.x(2:end), H(2:end), grid.x(2:end), flow.h(2:end,2)+flow.zb(2:end,2),...
%              grid.x(2:end),flow.h (2:end,2), grid.x(2:end), flow.zb(2:end,2),grid.x(2:end), v_h(2:end));
%         title(['time = ',num2str(itstep * run.dt),'s'])
%         legend('H','h+zb','h','zb','u²/2g','Location','north')
%         xlabel('x')
%         ylabel('[m]')
%         xlim([1 10]);
%         ylim([-0.1 1.2]);
%         hold off
%         pause(0.001)
        
        % Save Video
        drawnow
        currentFrame = getframe;
        Video(Frame) = getframe(gca);
        writeVideo(VidObj,currentFrame);    
        Frame = Frame +1;
        
        if itstep * run.dt == 0.05
        print(fig_Channeld,'-dpng',fullfile(dir,"time=0,05s"),'-r150');
        end
        
        if itstep * run.dt == 0.1
        print(fig_Channeld,'-dpng',fullfile(dir,"time=0,1s"),'-r150');
        end
               
        if itstep * run.dt == 1
        print(fig_Channeld,'-dpng',fullfile(dir,"time=1s"),'-r150');
        end
        
        if itstep * run.dt == 3
        print(fig_Channeld,'-dpng',fullfile(dir,"time=3s"),'-r150');
        end
        
        if itstep * run.dt == 5
        print(fig_Channeld,'-dpng',fullfile(dir,"time=5s"),'-r150');
        end
        
        if itstep * run.dt == 6
        print(fig_Channeld,'-dpng',fullfile(dir,"time=6s"),'-r150');
        end
        
        
     
        
    
    end
    
end


close(VidObj);
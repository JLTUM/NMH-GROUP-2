% **************************************************************************
% NMH Lab Summer Semester 2020, Assignment 9
%
% This code solves the 2D shallow water equations
%
% author: H. Zeng & L. Unglehrt
% July, 2020
% **************************************************************************
clear;
close all

%% Initialize simulation

global infilename 
infilename = "infile_2D_swe_damBreak_V1.mat";
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

% Create boundary conditions 
[ bconds ] = set_boundary_conditions();

%% Preallocation of variables

itdiag = 6;

U_channeld = [];
Q_channeld = [];

Frame = 1;
Video(run.ntst/itdiag) = struct('cdata',[],'colormap',[]);

VidObj = VideoWriter('kst30_ho1cm');  % Name of Video
open(VidObj)

disp(flow.kst(1,1))

%% Preallocation of Plots

 fig_Channeld = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
 Axes = gcf;
 xlim([1 10]);
 ylim([-0.1 1.2]);
 ax = gca;
 ax.NextPlot = 'replaceChildren';

%% Time integration

fprintf('start time integration\n')
for itstep = 1:run.ntst
    
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );
   
    if mod(itstep, itdiag) == 0

    % Diagonostic output
        
        % Energy and flow head
        u = flow.hu(:,2) ./ flow.h(:,2);                            % u velocity
        v = flow.hv(:,2) ./ flow.h(:,2);                            % v velocity
        U = sqrt (u.^2 + v.^2);                                     % Velocity Vector
        H = U.^2 / (2*constants.g) + flow.h(:,2) + flow.zb(:,2);    % energy head
        v_h = U.^2 / (2*constants.g);                               % velocity head
        c = sqrt(constants.g .* flow.h(:,2));                       % Wave celerity

        % CFL number
        fprintf('%d : CFL number:         %e\n', itstep, ...
            compute_CFL_number(constants, grid, run.dt, flow.h, flow.hu, flow.hv));
        
     % Plot results
    
        set(0, 'CurrentFigure', fig_Channeld)
%         plot(grid.x(2:end), H(2:end), grid.x(2:end), flow.h(2:end,2)+flow.zb(2:end,2),...
%              grid.x(2:end),flow.h (2:end,2), grid.x(2:end), flow.zb(2:end,2),grid.x(2:end), v_h(2:end));
%         title(['time = ',num2str(itstep * run.dt),'s'])
%         legend('H','h+zb','h','zb','u²/2g','Location','north')

        yyaxis left
        plot(grid.x(2:end), flow.h(2:end,2)+flow.zb(2:end,2),...
            grid.x(2:end), flow.zb(2:end,2),'k')
        xlim([1 10]);
        ylim([-0.1 1.2]);
%         yyaxis right
%         plot(grid.x(2:end), u(2:end), grid.x(2:end), c(2:end))
        
        title(['time = ',num2str(itstep * run.dt),'s',' kst = ', num2str(flow.kst(1,1))])
        legend('h+zb', 'Location','north')
    
%         s = surf(grid.x, grid.y, flow.h','EdgeColor', 'none','Facealpha','0.5');
%         view(15, 65)
%         xlim([1 10]);
%         ylim([-0.1 1.2]);

        xlabel('x')
        ylabel('[m]')

        hold off
        pause(0.001)
        
        % Save Video
        drawnow
        currentFrame = getframe;
        Video(Frame) = getframe(gca);
        writeVideo(VidObj,currentFrame);
        Frame = Frame +1;
        
    end
    
end

close(VidObj);

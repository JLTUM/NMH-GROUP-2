%**************************************************************************
% NMH Lab Summer Semester 2020, Assignment 7
%
% This code solves the 2D Shallow Water equations
%
% author: H. Zeng & L. Unglehrt
% June, 2020
%**************************************************************************

clear;
close all
set(0,'DefaultFigureVisible','off');
%% Initialize simulation
% read infile 
infilename = 'infile_2D_swe_test.mat';
fprintf('infilename is: %s\n', infilename)

% build structures 
[grid, run, constants, flow, bconds] = build_structs;
fprintf('struct built\n')

% fill some fields with data from input file
[grid, run, constants] = set_params(infilename);
fprintf('parameters set\n')

% Generate an equidistant grid 
[grid] = generate_grid(grid);    

fprintf('grid set\n')


% Set initial conditions 
run.t = 0;
[ flow ] = set_initial_condition( grid, flow );

% Create boundary conditions
bconds.bwest = {'WALL'};
bconds.beast = {'WALL'};
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

%% Time integration
for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );  
    
    if itstep == 1
        fig_Surf = figure('units','normalized','outerposition',[0 0 1 1]);
        fig_Quiver = figure('units','normalized','outerposition',[0 0 1 1]);
        surf(grid.x,grid.y,flow.h+flow.zb,'FaceAlpha',0.5)
        hold on
        surf(grid.x,grid.y,flow.zb,'FaceColor','b')
        xlabel('x','Fontsize',15)
        ylabel('y','Fontsize',15)
        zlabel('h','Fontsize',15)
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')
        set(gca,'YTickLabel',a,'fontsize',15,'FontWeight','bold')
        set(gca,'ZTickLabel',a,'fontsize',15,'FontWeight','bold')
        zlim([-1 2])
        title(['n=',num2str(itstep)])
        pause(0.05)
        hold off  
    end
    

%% Plot results

    set(0, 'CurrentFigure', fig_Surf) %Surf Definition
    surf(grid.x, grid.y, flow.h+flow.zb,'FaceAlpha',0.5)
    hold on 
    surf(grid.x, grid.y, flow.zb,'FaceColor','b')
    xlabel('x','Fontsize',15)
    ylabel('y','Fontsize',15)
    zlabel('h','Fontsize',15)
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')
    set(gca,'YTickLabel',a,'fontsize',15,'FontWeight','bold')
    set(gca,'ZTickLabel',a,'fontsize',15,'FontWeight','bold')
    zlim([-1 2])
    title(['n=',num2str(itstep)])
    pause(0.05)
    hold off

    set(0, 'CurrentFigure', fig_Quiver)%Quiver Definition
    quiver(grid.x,grid.y,flow.hu,flow.hv,'b')
    xlabel('x','Fontsize',15)
    ylabel('y','Fontsize',15)
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')
    set(gca,'YTickLabel',a,'fontsize',15,'FontWeight','bold')
    pause(0.05)


    if mod(itstep,10) == 0 || itstep == 1
        print(fig_Surf,'-dpng',sprintf("Plots_seven/Surf at n=%d.png", itstep),'-r150');
        print(fig_Quiver,'-dpng',sprintf("Plots_seven/Quiver at n=%d.png", itstep),'-r150');
    end


end

set(0,'DefaultFigureVisible','on');

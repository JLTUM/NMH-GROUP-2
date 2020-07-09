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
Fr = [];
diff = [];
% TODO: u_mean

<<<<<<< Updated upstream
%% Time integration
fprintf('start time integration\n')
for itstep = 1:run.ntst
=======
I = abs((flow.zb(2,2)-flow.zb(end))) / grid.xmax;
N_M = NWV_muster(Q,b,0,I,k_st);

%% Preallocation of Plots

    fig_Hy = figure;
        hold on
        title('H-y Diagram')
        y = 0:0.01:5;
        H = y + ( Q )^2 ./ (b * y ).^2 ./ 2 ./ 9.81;
        plot(H,y);
        H = y;
        y_c = (((Q/b)^2)/09.81)^(1/3);
        yline(y_c);
        plot(H,y);
        xlim([0 5])
        xlabel('H') 
        ylabel('y')
    fig_Channeld = figure('units','normalized','outerposition',[0 0 1 0.6]);
    
%% Time integration

for itstep = 1:40000 %run.ntst
>>>>>>> Stashed changes
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );

 %% Result treatments

    R_hyd(end+1) = (grid.ymax * nanmean(nanmean(flow.hu))) / (grid.ymax + (2 * nanmean(nanmean(flow.h))));
    v_st(end+1) = min(min(flow.kst)) * sqrt(abs(flow.I_S)) *  R_hyd(end)^(2/3);
    I_WSP(end+1) = abs((min(min(flow.h(2:end,:)))-max(max((flow.h(2:end,:))))))/grid.xmax;
    I_S(end+1) = abs((min(min(flow.zb))-max(max(flow.zb))))/grid.xmax;
    Fr(end+1) = mean(flow.hu(end-1,:) ./ sqrt( constants.g * flow.h(end-1,:).^3 ));
    % TODO: Froude number
    % TODO: u:mean
    
    % Channel Diagnosis
    %set(0, 'CurrentFigure', fig_D)
    
    
    
<<<<<<< Updated upstream
    figure(2)
    plot(1:itstep,I_S,'-b', 1:itstep,I_WSP,'-r', 1:itstep, v_st,'-g', 1:itstep, Fr,'-y')
    legend('I_S','I_{WSP}','v_{st}','Fr','Location','northwest')
    pause(0.01)
    
    % Water level 
    %set(0, 'CurrentFigure', fig_WSP)
    figure(3)
    plot(grid.x,flow.h(:,2),'b',grid.x,flow.hu(:,2),'g') %%plot(grid.x,flow.h(:,2)+flow.zb(:,2))
    legend('Waterdepth','Specific Discharge')
    title('Channel diagnosis')
    
    
    
if itstep == 1
         fig_Surf = figure('units','normalized','outerposition',[0 0 1 1]);
        fig_Quiver = figure('units','normalized','outerposition',[0 0 1 1]);
        surf(grid.x,grid.y,(flow.h+flow.zb)','FaceAlpha',0.5)
        hold on
        surf(grid.x,grid.y,flow.zb','FaceColor','b')
        xlabel('x','Fontsize',15)
        ylabel('y','Fontsize',15)
        zlabel('h','Fontsize',15)
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')
        set(gca,'YTickLabel',a,'fontsize',15,'FontWeight','bold')
        set(gca,'ZTickLabel',a,'fontsize',15,'FontWeight','bold')
        %zlim([-1 2])
        title(['n=',num2str(itstep)])
        pause(0.05)
        hold off  
    end
    

%% Plot results

   set(0, 'CurrentFigure', fig_Surf) %Surf Definition
    surf(grid.x, grid.y, (flow.h+flow.zb)','FaceAlpha',0.5)
    hold on 
    surf(grid.x, grid.y, flow.zb','FaceColor','b')
    xlabel('x','Fontsize',15)
    ylabel('y','Fontsize',15)
    zlabel('h','Fontsize',15)
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')
    set(gca,'YTickLabel',a,'fontsize',15,'FontWeight','bold')
    set(gca,'ZTickLabel',a,'fontsize',15,'FontWeight','bold')
    %zlim([-1 2])
    title(['n=',num2str(itstep)])
    pause(0.005)
    hold off
=======
%% Plot Results
  
if ~mod(itstep,500)
    % Plot Channel Diagnosis
    set(0, 'CurrentFigure', fig_Channeld)
%     subplot(131)
%     plot(grid.x(2:end), v_st(2:end), grid.x(2:end), u(2:end), grid.x(2:end), flow.hu(2:end,2));
%     legend('v_{st}','u','hu')
%     title('Channel Velocity')

    subplot(131)
    title('Channel Diagnosis')
    yyaxis left
    plot(1:itstep, Q_channeld)
    ylabel('[m³/s]')
    yyaxis right
    plot(1:itstep, U_channeld)
    ylabel('[m/s]')
    xlabel('timestep')
    legend('Q','U_b','Location','southwest')
    
    subplot(132)
    plot(grid.x(2:end), H(2:end), grid.x(2:end), flow.h(2:end,2)+flow.zb(2:end,2),...
         grid.x(2:end),flow.h (2:end,2), grid.x(2:end),flow.zb (2:end,2));
    title('Head and Waterlevel')
    legend('H','h+zb','h','zb','Location','southwest')
    xlabel('x')
    ylabel('[m]')    
    hold off
    
    subplot(133)
    title('Froud and specific Discharge')
    yyaxis left 
    plot(grid.x(2:end), Fr(2:end));
    ylabel('Froude [-]')
    yyaxis right
    plot(grid.x(2:end), flow.hu(2:end,2))
    ylabel('[m²/s]')
    xlabel('x')
    legend('Fr','hu','Location','northwest')

    sgtitle(['n= ',num2str(itstep)])

    % Plot Hy-Diagramm
%     set(0, 'CurrentFigure', fig_Hy)
%     hold on
%         scatter(H(2),flow.h(2,2),'b.')
%         scatter(H(976),flow.h(976,2),'g.')
%         scatter(H(1300),flow.h(1300,2),'r.')
%          legend(['x = ', num2str(grid.x(2))],['x = ', num2str(grid.x(976))],['x = ', num2str(grid.x(1300))])
%     
%          hold off
    
    pause(0.001)
    
end
>>>>>>> Stashed changes

%     set(0, 'CurrentFigure', fig_Quiver)%Quiver Definition
%     quiver(grid.y,grid.x,flow.hu,flow.hv,'b')
%     xlabel('x','Fontsize',15)
%     ylabel('y','Fontsize',15)
%     a = get(gca,'XTickLabel'); 
%     set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')
%     set(gca,'YTickLabel',a,'fontsize',15,'FontWeight','bold')
%     pause(0.05)

%     if mod(itstep,10) == 0 || itstep == 1
%         print(fig_Surf,'-dpng',sprintf("Plots_eight_Case1/Surf at n=%d.png", itstep),'-r150');
%         print(fig_Quiver,'-dpng',sprintf("Plots_eight_Case2/Quiver at n=%d.png", itstep),'-r150');
%     end    

end
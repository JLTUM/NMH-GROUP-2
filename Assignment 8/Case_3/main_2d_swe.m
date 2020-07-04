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
infilename = 'infile_2D_swe_channelFlow3.mat'; %% 1,2,3,4,5
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

%% Time integration
fprintf('start time integration\n')
for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );
    
%% PreProcessing 
%plot kst change over slope
% figure(999)
%     plot(grid.x,flow.zb(:,2),'m',grid.x(2:end-1),flow.kst(:,2),'r')
%     legend('bottom elevation','Strickler value (kst)','Location','northwest')
%     xlabel('x','Fontsize',15)
%     ylabel('y','Fontsize',15)
%     title('PreProcessing')
    
    

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
    
    
    
    figure(2)
    plot(1:itstep, v_st,'-g', 1:itstep, Fr,'-y')
    legend('v_{st}','Fr','Location','northwest')
    pause(0.01)
    
    figure(6)
    plot(1:itstep,I_S,'-b', 1:itstep,I_WSP,'-r')
    legend('I_S','I_{WSP}','Location','northwest')
    pause(0.01)
    
    % Water level 
    %set(0, 'CurrentFigure', fig_WSP)
    figure(3)
    plot(grid.x,flow.h(:,2),'b',grid.x,flow.hu(:,2),'r',grid.x,flow.zb(:,2),'g') %%plot(grid.x,flow.h(:,2)+flow.zb(:,2))
    legend('Waterdepth','Specific Discharge','bottom elevation')
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
if itst ==1
    fig_Quiver = figure('units','normalized','outerposition',[0 0 1 1]);
%     fig_Surf = figure('units','normalized','outerposition',[0 0 1 1]);
end
%% Quiver plot
% Grid
set(0, 'CurrentFigure', fig_Quiver)
figure(1)
quiver(grid.x(1:100),grid.y(1:100),flow.u,flow.v,'b');

%% Trace maximum of velo
[row, col] = find(ismember(flow.u, max(flow.u(:))));
figure(2)
surf(grid.x(1:100),grid.y(1:100),flow.u)
%plot(grid.x(col),grid.y(row),'+r')
hold on
%xlim([0,grid.x(grid.ntst)]);
%ylim([0,grid.y(grid.ntst)]);

%% Surf Plot
% fig_Surf = surf(grid.x(1:100),grid.y(1:100),flow.u);
% zlim([-1,1])


% if mod(itst,10) == 0
%     mkdir Plots_seven
%     print(fig_Quiver,'-dpng',sprintf("Plots_seven/Quiver_ntst=%d.png", itst),'-r150');
% 
%         figure(1)
%         subplot(2,1,1);
%         hold on
%         plot(grid.x(1:100),flow.u(51,:),"b")
%         plot(grid.x(1:100),flow_e.u(51,:),"g")
%         hold off
%         subplot(2,1,2);
%         hold on 
%         plot(grid.x(1:100),flow.v(51,:),"b")
%         plot(grid.x(1:100),flow_e.v(51,:),"g")
%         hold off
%         %end
%         temp = flow.u - flow_e.u;
%         average = mean(mean(temp));
%         error(itst,1) = mean(average);
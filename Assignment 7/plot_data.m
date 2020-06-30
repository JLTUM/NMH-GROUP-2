if itst ==1
%     fig_Quiver = figure('units','normalized','outerposition',[0 0 1 1]);
    fig_Surf = figure('units','normalized','outerposition',[0 0 1 1]);
end
%% Quiver plot
% % Grid
% X = ndgrid(grid.x(1:size(flow.u,2)));
% Y = ndgrid(grid.y(1:size(flow.u,2)));
% set(0, 'CurrentFigure', fig_Quiver)
% quiver(X,Y,flow.u,flow.v,'b');
% xlabel('x','Fontsize',15)
% ylabel('y','Fontsize',15)
% title(['Velocity field at n = ',num2str(itst)])

%% Trace maximum of velo
% max_value = max(max(flow.u));
% loc = find(abs(flow.u-max_value)<5e-4);
% [row,col] = ind2sub(size(flow.u),loc);
% plot(grid.x(col),grid.y(row),'+r')
% hold on
% hold off
% hold on
% xlim([0,10]);
% ylim([0,10]);

% Surf Plot
set(0, 'CurrentFigure', fig_Surf)
surf(grid.x(1:size(flow.u,2)),grid.y(1:size(flow.u,2)),flow.u);
zlim([-1,1])


% if mod(itst,100) == 0
%     mkdir Plots_seven
%     print(fig_Quiver,'-dpng',sprintf("Plots_seven/Quiver_ntst=%d.png", itst),'-r150');
% end
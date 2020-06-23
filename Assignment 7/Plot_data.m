if itst ==1
    fig_Quiver = figure('units','normalized','outerposition',[0 0 1 1]);
    fig_Surf = figure('units','normalized','outerposition',[0 0 1 1]);
end


set(0, 'CurrentFigure', fig_Surf)
surf(grid.x(1:grid.nx),grid.y(1:grid.ny),flow.u);
title("Surface Velocity Plot of Runge-Kutta Approximation")
xlabel("x",'FontName','Times','fontsize',18)
ylabel("y",'FontName','Times','fontsize',18)
zlabel("u",'FontName','Times','fontsize',18)
pause(0.005)


%% Quiver plot
%Grid
set(0, 'CurrentFigure', fig_Quiver)
X = ndgrid(grid.x(1:grid.nx));
Y = ndgrid(grid.y(1:grid.ny));
quiver(X,Y,flow.u,flow.v,'b');
xlabel("x",'FontName','Times','fontsize',18)
ylabel("y",'FontName','Times','fontsize',18)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'FontName','Times','fontsize',18)

%% Print
if mod(itst,100) == 0
    mkdir Plots_seven
    title(['Velocity Field of Runge-Kutta Timestep: ',num2str(itst)],'FontName','Times','fontsize',18)
    print(fig_Quiver,'-dpng',sprintf("Plots_seven/Quiver_ntst=%d.png", itst),'-r150');
end

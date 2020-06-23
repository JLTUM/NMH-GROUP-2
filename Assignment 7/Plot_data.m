if itst ==1
    fig_Quiver = figure('units','normalized','outerposition',[0 0 1 1]);
    fig_Surf = figure('units','normalized','outerposition',[0 0 1 1]);
end

set(0, 'CurrentFigure', fig_Surf)
surf(grid.x(1:grid.nx),grid.y(1:grid.ny),flow.u);
title("Surface Velocity Plot of Runge-Kutta Approximation")
xlabel("x")
ylabel("y")
zlabel("u")
pause(0.005)


%% Quiver plot
%Grid
set(0, 'CurrentFigure', fig_Quiver)
X = ndgrid(grid.x(1:grid.nx));
Y = ndgrid(grid.y(1:grid.ny));
quiver(X,Y,flow.u,flow.v,'b');
title("Velocity Field of Runge-Kutta Approximation")
xlabel("x")
ylabel("y")


%% Print
if mod(itst,100) == 0
    mkdir Plots_seven
    print(fig_Quiver,'-dpng',sprintf("Plots_seven/Quiver_ntst=%d.png", itst),'-r150');
end

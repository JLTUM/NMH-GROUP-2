% if itst ==1
%     fig_Quiver = figure('units','normalized','outerposition',[0 0 1 1]);
%     fig_Surf = figure('units','normalized','outerposition',[0 0 1 1]);
% end

% figure(1)
% quiver(grid.x(1:100),grid.y(1:100),flow.u,flow.v,'b');
% 
% [row, col] = find(ismember(flow.u, max(flow.u(:))));
%plot(grid.x(col),grid.y(row),'+r')
figure(1)
surf(grid.x(1:30),grid.y(1:30),flow.u);
pause(0.005)

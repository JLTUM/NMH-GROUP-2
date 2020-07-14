grid.x = 1:0.1:100;
grid.y = 1:1:4;


flow.h = zeros( length(grid.x), length(grid.y) );
    flow.h(:,:) = 0.001;
    flow.h(1:51,:) = 1;
    
    
    plot(grid.x, flow.h)
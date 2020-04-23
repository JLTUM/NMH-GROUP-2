clear all
%these are Julians comments 
for exp = 1:1:4

    % Number of grid points
    n = 10^exp;

    % Grid spacing
    h = 2.0*pi / (n-1);

    % Calculation of... 
    for i = 1 : n

      % ...coordinates
      x(i) = 0.0 + (i-1) * h;

      % ...analytical values of function and derivative
      f(i) = sin(x(i));
      dfe(i) = cos( x(i) );
    end

    for i = 1 : n

      % Periodic boundary conditions... 

      % ...for the first point
      if ( i == 1 )

        fp = f(i+1);
        fm = f(n-1);

      end

      % ...for the last point
      if ( i == n )

        fp = f(2);
        fm = f(n-1);

      end


      % If point is not at the boundaries 
      if ( ( i ~= 1 ) && ( i ~= n ) )
        fp = f(i+1);
        fm = f(i-1);
      end

      fi = f(i);

    % First derivative: Upwind
        dfn_U(i) = ( fi - fm ) / ( 1.0*h );
      
    % First derivative: Downwind
        dfn_D(i) = ( fp - fi ) / ( 1.0*h );
        
    % First derivative: Central Limit
        dfn_C(i) = ( fp - fm) / ( 2.0*h );

    % Error: Upwind
        er_U(i) = abs( ( dfe(i) - dfn_U(i) ) / dfe(i) );
    
    % Error: Upwind
        er_D(i) = abs( ( dfe(i) - dfn_D(i) ) / dfe(i) );
      
    % Error: Upwind
        er_C(i) = abs( ( dfe(i) - dfn_C(i) ) / dfe(i) );
    
      
    end
    
    figure
    % Plotting of analytical solution and numerical approximation
    plot(x, dfn_U, x, dfn_U, x, dfn_C, x, dfe)
    legend('Upwind','Downwind','Central Limit','Exact','Location','SouthEast')
    set(gca,'FontSize',14); 
    title(n);

    % Storing grid spacing and error in matrix 'error'
    error(exp,1) = h;
    error(exp,2) = er_U(n/5);
    error(exp,3) = er_D(n/5);
    error(exp,4) = er_C(n/5);
    

end


% Plotting of error over grid spacing in normal scale
loglog(error(:,1),error(:,2),error(:,1),error(:,3),error(:,1),error(:,4));
%set(gca,'FontSize',14); 
title('Error plot');
xlabel('Grid spacing h');
ylabel('Relative error');
legend('Error Upwind','Error Downwind','Error Central Limit','Location','SouthEast')

figure

% Changing axes to logarithmic scale 
loglog(error(:,1),error(:,2),'-ro');
%set(gca,'FontSize',14); 
title('Error plot');
xlabel('Grid spacing h');
ylabel('Relative error');
legend('Error','Location','SouthEast')

figure

% Adding the linear function for comparison 
loglog(error(:,1),error(:,2),'-ro',x,x);
%set(gca,'FontSize',14); 
title('Error plot');
xlabel('Grid spacing h');
ylabel('Relative error');
legend('Error','x','Location','SouthEast')

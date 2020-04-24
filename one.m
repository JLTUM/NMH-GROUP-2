%% Header

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 1
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

clear all
close all

%% grid and analytical solution 

for exp = 1:1:4

    % Number of grid points
    n = 10^exp

    % Grid spacing
    h = 2.0*pi / (n-1);

    % Calculation of... 
    for i = 1 : n

      % ...coordinates
      x(i) = 0.0 + (i-1) * h;

      % ...analytical values of function and derivative
      f(i) = sin(x(i)); % analytical
      dfe(i) = cos( x(i) ); % analytical derivative
    end
    
%% Discretization for numerical solution

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

%% function values 

    % First derivative: Upwind
        dfn_U(i) = ( fi - fm ) / ( 1.0*h );
      
    % First derivative: Downwind
        dfn_D(i) = ( fp - fi ) / ( 1.0*h );
        
    % First derivative: Central
        dfn_C(i) = ( fp - fm) / ( 2.0*h );

    % Error: Upwind
        er_U(i) = abs( ( dfe(i) - dfn_U(i) ) / dfe(i) );
    
    % Error: Downwind
        er_D(i) = abs( ( dfe(i) - dfn_D(i) ) / dfe(i) );
      
    % Error: Central
        er_C(i) = abs( ( dfe(i) - dfn_C(i) ) / dfe(i) );
    
      
    end

    % Plotting of analytical solution and numerical approximation of first
    % derivatives
    subplot(2,2,exp)
    plot(x, dfn_U, x, dfn_U, x, dfn_C, x, dfe)
    legend('Upwind','Downwind','Central','Exact')
    set(gca,'FontSize',14); 
    title(n);

    % Storing grid spacing and error in matrix 'error'
    error(exp,1) = h;
    error(exp,2) = er_U(n/5);
    error(exp,3) = er_D(n/5);
    error(exp,4) = er_C(n/5);
    
end

figure
% Plotting of error over grid spacing in normal scale
loglog(error(:,1),error(:,2),'-x',error(:,1),error(:,3),'-x',error(:,1),error(:,4),'-x', x, x,'-k');
title('Error plot');
xlabel('Grid spacing h');
ylabel('Relative error');
legend('Error Upwind','Error Downwind','Error Central','Location','SouthEast')

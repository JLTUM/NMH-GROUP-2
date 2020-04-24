%% Header

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 1
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

clear all
close all

%% grid and analytical solution 

grid = [10,20,50,100,1000,10000];
% grid = linspace(10,10000,1000);
for j = 1:length(grid)

    % Number of grid points
    n = grid(j);

    % Grid spacing
    h = 2.0*pi / (n-1);

    % Calculation of... 
    for i = 1 : n

      % ...coordinates
      x(i) = 0.0 + (i-1) * h;

      % ...analytical values of function and derivative
      f(i) = sin(2*x(i)); % analytical
      dfe(i) = cos( 2*x(i) ) * 2; % analytical first derivative
      dfe2(i) = -sin( 2*x(i) ) * 4; % analytical second derivative
      
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
        dfn_C(i) = ( fp - fm ) / ( 2.0*h );
        
    % Second Derivative: Central 
        dfn_C2(i) = ( fp - 2*fi + fm ) / (h^2);

    % Error: Upwind
        er_U(i) = abs( ( dfe(i) - dfn_U(i) ) / dfe(i) );
    
    % Error: Downwind
        er_D(i) = abs( ( dfe(i) - dfn_D(i) ) / dfe(i) );
      
    % Error: Central first derivative
        er_C(i) = abs( ( dfe(i) - dfn_C(i) ) / dfe(i) );
        
    % Error: Central second derivative
        er_C2(i) = abs ( ( dfe2(i) - dfn_C2(i) ) / dfe(i) );
        
        
    end  

%% Interpolation of Points

    xq = [0+h/2:h:n];
    dfn_U_interp = interp1(x,dfn_U,xq);
    dfn_D_interp = interp1(x,dfn_D,xq);
    dfn_C_interp = interp1(x,dfn_C,xq);
    dfn_C2_interp = interp1(x,dfn_C2,xq);
    
    %% Plots
    
    if j == 1 
        fig_first_derivative = figure;
        fig_second_derivative = figure;
    end
    
    % plot first derivative 
    set(0,'CurrentFigure',fig_first_derivative)
    subplot(3,2,j)
    plot(x, dfn_U, 'b.:', x, dfn_D, 'r.:', x, dfn_C, 'g.:', x, dfe, '--k') % numerical values
    hold on 
    plot(xq, dfn_U_interp, 'b.', xq, dfn_D_interp, 'r.', xq, dfn_C_interp, 'g.') % interpolated values
    legend('Upwind','Downwind','Central', 'Exact')
    set(gca,'FontSize',14); 
    title(n);

    % plot second derivative
    set(0,'CurrentFigure',fig_second_derivative)
    subplot(3,2,j)
    plot(x, dfn_C2, 'g.:', x, dfe2, '--k')
    hold on
    plot(xq, dfn_C2_interp, 'g.')
    legend('Central', 'Exact')
    set(gca,'FontSize',14); 
    title(n);
    
    % Storing grid spacing and error in matrix 'error'
    error(j,1) = h;
    error(j,2) = er_U(n/5);
    error(j,3) = er_D(n/5);
    error(j,4) = er_C(n/5);
    error(j,5) = er_C2(n/5);
    
end
% Plotting of error over grid spacing in normal scale
figure
    loglog(error(:,1),error(:,2),'-',error(:,1),error(:,3),'-',error(:,1),error(:,4),'-', x, x,'-k');
    title('Error plot first derivative');
    xlabel('Grid spacing h');
    ylabel('Relative error');
    legend('Error Upwind','Error Downwind','Error Central','Location','SouthEast')

% Plotting of error over grid spacing in normal scale
figure
    loglog(error(:,1),error(:,5),'-', x, x,'-k');
    title('Error plot second derivative');
    xlabel('Grid spacing h');
    ylabel('Relative error');
    legend('Error Central','Location','SouthEast')

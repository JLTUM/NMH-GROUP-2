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
    n = 10^exp;

    % Grid spacing
    h = 2.0*pi / (n-1);

    % Calculation of... 
    for i = 1 : n

      % ...coordinates
      x(i) = 0.0 + (i-1) * h;

      % ...analytical values of function and derivative
      f(i) = sin( 2*x(i) ); % analytical
      dfe(i) = 2*cos( 2*x(i) ); % analytical derivative
      ddfe(i) = -4*sin ( 2*x(i) ); %second derivative
      
      %changed these to 2pi instead of pi because we should model based on
      %group number
      
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
    
        
    % Second Derivative Central:
        ddfn_C(i) = (fp - 2*fi + fm)/(h^2);
    % Error:
        er_DC(i) = abs( ( ddfe(i) - ddfn_(i) ) / ddfe(i) );
    
      
    end  
    % -----First Derivative -----
    % Plotting of analytical solution and numerical approximation of first
    % derivatives
    subplot(4,4,exp)
    plot(x, dfn_U, x, dfn_D, x, dfn_C, x, dfe)
    legend('Upwind','Downwind','Central','Exact')
    set(gca,'FontSize',14); 
    title("1st" + n);

    % Storing grid spacing and error in matrix 'error'
    error(exp,1) = h;
    error(exp,2) = er_U(n/5);
    error(exp,3) = er_D(n/5);
    error(exp,4) = er_C(n/5);

    % -----Second Derivative -----
    % Plotting of analytical solution and numerical approximation of first
    % derivatives
    subplot(4,4,exp+4)
    plot(x, ddfn_C, x, ddfe)
    legend('Approx','Exact')
    set(gca,'FontSize',14); 
    title("2nd" + n);
    %issues with this method of representation think about making separate
    %loops for 2nd
    % Storing grid spacing and error in matrix 'error'
    error(exp,5) = er_DC(n/5);
    
end

figure
% Plotting of error over grid spacing in normal scale
loglog(error(:,1),error(:,2),'-x',error(:,1),error(:,3),'-x',error(:,1),error(:,4),'-x',error(:,1),error(:,5),'-x',x, x,'-k');
title('Error plot');
xlabel('Grid spacing h');
ylabel('Relative error');
legend('Error Upwind','Error Downwind','Error Central',"Error 2nd Derivative",'Location','SouthEast')

%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 1
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

% ----------------------------------------------------------------------- %

clear all
close all

%% grid and analytical solution 

grid = [10,20,50,100,1000,10000]; % amount of grid points 
 
for j = 1:length(grid)

    % Number of grid points
    n = grid(j);

    % Grid spacing
    h = 2.0 * pi / (n-1);

    % Calculation of... 
    for i = 1 : n

        % ...coordinates
        x(i) = 0.0 + (i-1) * h;

        % ...analytical values of function and derivative
        f(i) = sin(2*x(i)); % analytical function
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
    
    xq = [0+h/2:h:h*(n-1)]; % x-values on which function values get interpolated: points in the middle of two grid points
    dfn_U_interp = interp1(x,dfn_U,xq); % Interpolated Points for Upwind
    dfn_D_interp = interp1(x,dfn_D,xq); % Interpolated Points for Downwind
    dfn_C_interp = interp1(x,dfn_C,xq); % Interpolated Points for Central
    dfn_C2_interp = interp1(x,dfn_C2,xq); % Interpolated Points for Central (second der.)
    
    for k = 1:n-1
    dfe_int(k) = cos( 2*xq(k) ) * 2; % analytical first derivative
    dfe2_int(k) = -sin( 2*xq(k) ) * 4; % analytical second derivative
  
        
        %% function values 

    % Error: Int_Upwind
        int_er_U(k) = abs( ( dfe_int(k) - dfn_U_interp(k) ) / dfe_int(k) );
    
    % Error: Int_Downwind
        int_er_D(k) = abs( ( dfe_int(k) - dfn_D_interp(k) ) / dfe_int(k) );
      
    % Error: Central first derivative
        int_er_C(k) = abs( ( dfe_int(k) - dfn_C_interp(k) ) / dfe_int(k) );
        
    % Error: Central second derivative
        int_er_C2(k) = abs ( ( dfe2_int(k) - dfn_C2_interp(k) ) / dfe_int(k) ); 
        
    end
    
    
    
    
%% Plots
    
    if j == 1 
        fig_first_derivative = figure;
        fig_second_derivative = figure;
    end
    
%     %plot first derivative 
%     set(0,'CurrentFigure',fig_first_derivative)
%     subplot(3,2,j)
%     %plot(x, dfn_U, '-', x, dfn_D, '-', x, dfn_C, '-', x, dfe, '-k') % values of numerical difference scheme
%     hold on 
%     plot(xq, dfn_U_interp, '-', xq, dfn_D_interp, '-', xq, dfn_C_interp, '-',x, dfe, '-k') % interpolated values
%     legend('Upwind','Downwind','Central', 'Exact')
%     set(gca,'FontSize',14); 
%     title(n);
%     hold off

    % plot second derivative
    set(0,'CurrentFigure',fig_second_derivative)
    subplot(1,1,j)
    plot(x, dfn_C2, '-x', x, dfe2, '-k') % values of numerical difference scheme
    hold on
    plot(xq, dfn_C2_interp, '-o') % interpolated values
    legend('Central', 'Exact', "Interpolated")
    set(gca,'FontSize',14); 
    title(n);
    
    % Storing grid spacing and relative error in matrix
    error(j,1) = h;
    error(j,2) = er_U(n/5);
    error(j,3) = er_D(n/5);
    error(j,4) = er_C(n/5);
    error(j,5) = er_C2(n/5);
    error(j,6) = int_er_U(n/5);
    error(j,7) = int_er_D(n/5);
    error(j,8) = int_er_C(n/5);
    error(j,9) = int_er_C2(n/5);
    error(j,10)=n;
    
    
    
end

% Plotting of error over grid spacing in log scale
figure
    loglog(error(:,1),error(:,2),'-x',error(:,1),error(:,3),'-x',error(:,1),error(:,4),'-x'...
        ,error(:,1),error(:,6),'-o',error(:,1),error(:,7),'-o',error(:,1),error(:,8),'-o', x, x,'-k');
    title('Error plot first derivative');
    xlabel('Grid spacing h');
    ylabel('Relative error');
    legend('Error Upwind','Error Downwind','Error Central','Error Upwind int.','Error Downwind int.','Error Central int.','Location','SouthEast')

% Plotting of error over grid spacing in log scale

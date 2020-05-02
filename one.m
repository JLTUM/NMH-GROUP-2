%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 1
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

% ----------------------------------------------------------------------- %

clear all
close all

%% grid and analytical solution 

%grid = [10,12,13,14,15,16,17,18,20,30,50,100,1000,10000];
grid = [10,20,30,50,70,100,1000,10000];
% grid = [10:5:1000];
%grid = [2,4,6,8,10,12,14,16,18,20,30,40,60,100,1000];
% amount of grid points 
 
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
        f(i) = sin( 2*x(i) ); % analytical function
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

%% Numerical shemes

    % First derivative: Upwind
    dfn_U(i) = ( fi - fm ) / ( 1.0*h );
      
    % First derivative: Downwind
    dfn_D(i) = ( fp - fi ) / ( 1.0*h );
        
    % First derivative: Central
    dfn_C(i) = ( fp - fm ) / ( 2.0*h );
        
    % Second Derivative: Central 
    dfn_C2(i) = ( fp - 2*fi + fm ) / (h^2);
      
    end
    
%% Error of numerical shemes
    
    % Location of Error 
    x_error = 2*pi/5;
    
    % Value at Location of Error
    val_error_U = interp1(x,dfn_U,x_error);     % Value at linearily interpolated locatin of error: Upwind
    val_error_D = interp1(x,dfn_D,x_error);     % Value at linearily interpolated locatin of error: Downwind       
    val_error_C = interp1(x,dfn_C,x_error);     % Value at linearily interpolated locatin of error: Central
    val_error_C2 = interp1(x,dfn_C2,x_error);   % Value at linearily interpolated locatin of error: Central second derivative
    val_error_dfe = cos( 2*x_error ) * 2;       % Value at location of error first derivative
    val_error_dfe2 = -sin( 2*x_error ) * 4;     % Value at location of error second derivative
    
    % Error 
    er_U = abs( ( val_error_dfe - val_error_U ) / val_error_dfe ) ;
    er_D = abs( ( val_error_dfe - val_error_D ) / val_error_dfe ) ;
    er_C = abs( ( val_error_dfe - val_error_C ) / val_error_dfe ) ;
    er_C2 = abs( ( val_error_dfe2 - val_error_C2 ) / val_error_dfe2 ) ;
    
%% Interpolation of Points
    
    % Location of interpolated points
    xq = [0+h/2:h:h*(n-1)];               % x-values on which function values get interpolated

    % Value at interpolated points
    dfn_U_interp = interp1(x,dfn_U,xq);   % Interpolated Points for Upwind
    dfn_D_interp = interp1(x,dfn_D,xq);   % Interpolated Points for Downwind
    dfn_C_interp = interp1(x,dfn_C,xq);   % Interpolated Points for Central
    dfn_C2_interp = interp1(x,dfn_C2,xq); % Interpolated Points for Central (second der.)
    
%% Error of interpolated points

    % Value at Location of Error
    val_error_U_interp = interp1(xq,dfn_U_interp,x_error);   % Value at linearily interpolated locatin of error: Upwind
    val_error_D_interp = interp1(xq,dfn_D_interp,x_error);   % Value at linearily interpolated locatin of error: Downwind       
    val_error_C_interp = interp1(xq,dfn_C_interp,x_error);   % Value at linearily interpolated locatin of error: Central
    val_error_C2_interp = interp1(xq,dfn_C2_interp,x_error); % Value at linearily interpolated locatin of error: Central second derivative
 
    % Error
    er_U_interp = abs( ( val_error_dfe - val_error_U_interp ) / val_error_dfe ) ;
    er_D_interp = abs( ( val_error_dfe - val_error_D_interp ) / val_error_dfe ) ;
    er_C_interp = abs( ( val_error_dfe - val_error_C_interp ) / val_error_dfe ) ;
    er_C2_interp = abs( ( val_error_dfe2 - val_error_C2_interp ) / val_error_dfe2 ) ;
    
%% Error matrix

    % Error matrix
    error(j,1) = h;            % grid spacing
    error(j,2) = er_U;                  
    error(j,3) = er_D;             
    error(j,4) = er_C;             
    error(j,5) = er_U_interp;   
    error(j,6) = er_D_interp;   
    error(j,7) = er_C_interp;   
    error(j,8) = er_C2;            
    error(j,9) = er_C2_interp;  
        
%% Plots
    
    if j == 1 
        fig_first_derivative = figure;
        fig_second_derivative = figure;
    end
    
%     plot first derivative 
    set(0,'CurrentFigure',fig_first_derivative)
    subplot(4,2,j)
    plot(x, dfn_U, 'x-', x, dfn_D, 'x-', x, dfn_C, 'x-', x, dfe, '-k') % values of numerical difference scheme
    hold on 
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(xq, dfn_U_interp, '.--', xq, dfn_D_interp, '.--', xq, dfn_C_interp, '.--') % interpolated values
    legend('Upwind','Downwind','Central', 'Exact','Autoupdate','off')
    set(gca,'FontSize',14); 
    title(n);
    xline(x_error);
    hold off

%     plot second derivative
    set(0,'CurrentFigure',fig_second_derivative)
    subplot(4,2,j)
    ax = gca;
    ax.ColorOrderIndex = 4;
    hold on
    plot(x, dfn_C2, '-x', x, dfe2, '-k') % values of numerical difference scheme
    ax.ColorOrderIndex = 4;
    plot(xq, dfn_C2_interp, '.--') % interpolated values
    legend('Central', 'Exact','Autoupdate','off')
    set(gca,'FontSize',14); 
    xline(x_error);
    title(n);
    
end
    
% Plotting of error over grid spacing in log scale
figure
    loglog(error(:,1),error(:,2),'-',error(:,1),error(:,3),'-',...
        error(:,1),error(:,4),'-',error(:,1),error(:,8),'-')
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    loglog(error(:,1),error(:,5),'--',error(:,1),error(:,6),'--',...
        error(:,1),error(:,7),'--',error(:,1),error(:,9),'--')
    title('Error plot first derivative');
    xlabel('Grid spacing h');
    ylabel('Relative error');
    legend('Upwind','Downwind','Central','Central Second',...
        'Upwind Int','Downwind Int','Central Int','Central Second Int',...
        'Location','SouthEast')
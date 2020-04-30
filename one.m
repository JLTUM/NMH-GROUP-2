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

%% function values and errors

    % First derivative: Upwind
    dfn_U(i) = ( fi - fm ) / ( 1.0*h );
      
    % First derivative: Downwind
    dfn_D(i) = ( fp - fi ) / ( 1.0*h );
        
    % First derivative: Central
    dfn_C(i) = ( fp - fm ) / ( 2.0*h );
        
    % Second Derivative: Central 
    dfn_C2(i) = ( fp - 2*fi + fm ) / (h^2);

    % Error: First derivative Upwind
    er_U(i) = abs( ( dfe(i) - dfn_U(i) ) / dfe(i) );
    
    % Error: First derivative Downwind
    er_D(i) = abs( ( dfe(i) - dfn_D(i) ) / dfe(i) );
      
    % Error: First derivative Central
    er_C(i) = abs( ( dfe(i) - dfn_C(i) ) / dfe(i) );
        
    % Error: Second derivative Central
    er_C2(i) = abs ( ( dfe2(i) - dfn_C2(i) ) / dfe(i) );
      
    end
    
%% Interpolation of Points and error on interpolated points
    
    xq = [0:h/2:h*(n-1)]; % x-values on which function values get interpolated

    dfn_U_interp = interp1(x,dfn_U,xq); % Interpolated Points for Upwind
    dfn_D_interp = interp1(x,dfn_D,xq); % Interpolated Points for Downwind
    dfn_C_interp = interp1(x,dfn_C,xq); % Interpolated Points for Central
    dfn_C2_interp = interp1(x,dfn_C2,xq); % Interpolated Points for Central (second der.)
    
    dfe_interp = cos( 2 .* xq ) * 2; % analytical first derivative
    dfe2_interp = -sin( 2 .* xq ) * 4; % analytical second derivative
    
%% Error with interpolated points
    
    % Error: Int_Upwind
    int_er_U = abs( ( dfe_interp - dfn_U_interp ) ./ dfe_interp );
    
    % Error: Int_Downwind
    int_er_D = abs( ( dfe_interp - dfn_D_interp ) ./ dfe_interp );
      
    % Error: Central first derivative
    int_er_C = abs( ( dfe_interp - dfn_C_interp ) ./ dfe_interp );
        
    % Error: Central second derivative
    int_er_C2 = abs ( ( dfe2_interp - dfn_C2_interp ) ./ dfe_interp ); 
    
%% Plots
    
    if j == 1 
        fig_first_derivative = figure;
        fig_second_derivative = figure;
    end
    
%     plot first derivative 
    set(0,'CurrentFigure',fig_first_derivative)
    subplot(4,2,j)
    plot(x, dfn_U, 'o-', x, dfn_D, 'o-', x, dfn_C, 'o-', x, dfe, '--k') % values of numerical difference scheme
    hold on 
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(xq, dfn_U_interp, '.-', xq, dfn_D_interp, '.-', xq, dfn_C_interp, '.-') % interpolated values
    legend('Upwind','Downwind','Central', 'Exact')
    set(gca,'FontSize',14); 
    title(n);
    hold off

%     plot second derivative
    set(0,'CurrentFigure',fig_second_derivative)
    subplot(4,2,j)
    ax = gca;
    ax.ColorOrderIndex = 4;
    hold on
    plot(x, dfn_C2, '-o', x, dfe2, '--k') % values of numerical difference scheme
    ax.ColorOrderIndex = 4;
    plot(xq, dfn_C2_interp, '.-') % interpolated values
    legend('Central', 'Exact')
    set(gca,'FontSize',14); 
    title(n);
    
    % Storing grid spacing and relative error in matrix
%         error(j,1) = h;
%         error(j,2) = mean(er_U);
%         error(j,3) = mean(er_D);
%         error(j,4) = mean(er_C);
%         error(j,5) = mean(er_C2);
%         error(j,6) = mean(int_er_U);
%         error(j,7) = mean(int_er_D);
%         error(j,8) = mean(int_er_C);
%         error(j,9) = mean(int_er_C2);

    error(j,1) = h;
%     error(j,1) = n;
    error(j,2) = er_U(n/5); % value of 2*pi/5 for x
    error(j,3) = er_D(n/5); % value of 2*pi/5 for x
    error(j,4) = er_C(n/5); % value of 2*pi/5 for x
    error(j,5) = int_er_U((2*n/5)-1); % value of 2*pi/5 for xq
    error(j,6) = int_er_D((2*n/5)-1); % value of 2*pi/5 for xq
    error(j,7) = int_er_C((2*n/5)-1); % value of 2*pi/5 for xq
    error(j,8) = er_C2(n/5); % value of 2*pi/5 for x
    error(j,9) = int_er_C2((2*n/5)-1); % value of 2*pi/5 for xq
    
    end
    
% Plotting of error over grid spacing in log scale
figure
    loglog(error(:,1),error(:,2),'-',error(:,1),error(:,3),'-',error(:,1),error(:,4),'-')
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    loglog(error(:,1),error(:,5),'--',error(:,1),error(:,6),'--',error(:,1),error(:,7),'--')
    title('Error plot first derivative');
    xlabel('Grid spacing h');
    ylabel('Relative error');
    legend('Upwind','Downwind','Central',...
        'Upwind Int','Downwind Int','Central Int',...
        'Location','SouthEast')

% Plotting of error over grid spacing in log scale
figure
    ax = gca;
    ax.ColorOrderIndex = 4;
    hold on
    loglog(error(:,1),error(:,8),'-')
    hold on
    ax.ColorOrderIndex = 4;
    loglog(error(:,1),error(:,9),'--')
    ax.XScale = 'Log';
    ax.YScale = 'Log';
    title('Error plot second derivative');
    xlabel('Grid spacing h');
    ylabel('Relative error');
    legend('Central 2nd','Central 2nd Int','Location','SouthEast')
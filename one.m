%% ---------------------------- Header --------------------------------- %%

%%% Numerische Methoden der Hydromechanik
%%% Assignment: 1
%%% Group: 2
%%% Members: Nick Pfeiffer, Andreas Mirlach, Julian Lenz, Faro Schäfer

% ----------------------------------------------------------------------- %

clear all
close all
addpath('Plots_one')

%% grid and analytical solution 

grid = [10,20,30,50,100,1000,10000];

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
    val_error_U_interp = interp1(xq,dfn_U_interp,x_error);   % Value at linearely interpolated locatin of error: Upwind
    val_error_D_interp = interp1(xq,dfn_D_interp,x_error);   % Value at linearely interpolated locatin of error: Downwind       
    val_error_C_interp = interp1(xq,dfn_C_interp,x_error);   % Value at linearely interpolated locatin of error: Central
    val_error_C2_interp = interp1(xq,dfn_C2_interp,x_error); % Value at linearely interpolated locatin of error: Central second derivative
 
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
        % set X-tick Lables
            xticks = [0, 2*pi/5, pi/2, pi, 3*pi/2, 2*pi];
            xticklabels = {'\fontsize{11} 0','\color[rgb]{0.5,0.5,0.5}\fontsize{11}2\pi/5','\fontsize{11} \pi/2','\fontsize{11} \pi','\fontsize{11} 3\pi/2','\fontsize{11} 2\pi'};
        fig_first_derivative = figure('units','normalized','outerposition',[0 0 1 1]);
        fig_second_derivative = figure('units','normalized','outerposition',[0 0 1 1]);
        fig_first_derivative_int = figure('units','normalized','outerposition',[0 0 1 1]);

    end

    if (6<j)            % no need to plot high n (look identical to n = 1000) 
        continue
    end
    
%     plot first derivative 
    set(0,'CurrentFigure',fig_first_derivative)
    subplot(3,2,j)
    plot(x, dfn_U, '-', x, dfn_D, '-', x, dfn_C, '-', x, dfe, '-k') % values of numerical difference scheme
    hold on 
    ax = gca;
    ax.ColorOrderIndex = 1;
    legend('Upwind','Downwind', 'Central','Exact','Autoupdate','off','Location','best')
    set(gca,'FontSize',14); 
    xlim([0,2*pi])
    title(['First Derivative n=',num2str(n)]);
    set(gca,'XTick',xticks)
    ax.XTickLabel = xticklabels;
    hold off

%     plot second derivative
    set(0,'CurrentFigure',fig_second_derivative)
    subplot(3,2,j)
    ax = gca;
    ax.ColorOrderIndex = 5;
    hold on
    plot(x, dfe2, '-k',x, dfn_C2, '-') % values of numerical difference scheme
    ax.ColorOrderIndex = 5;
    plot(xq, dfn_C2_interp, '.--') % interpolated values
    legend('Exact','Central','Central Int','Autoupdate','off','Location','best')
    set(gca,'FontSize',14); 
    xlim([0,2*pi])
    title(['Second Derivative n=',num2str(n)]);
    set(gca,'XTick',xticks)
    ax.XTickLabel = xticklabels;  
    hold off
    
    if (j>2)      
        continue
    end
    
% plot first derivative interpolated
    set(0,'CurrentFigure',fig_first_derivative_int)

% Upwind
    subplot(3,2,j)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    xlim([0,2*pi])
    title(['n=',num2str(n)]);
    set(gca,'XTick',xticks)
    ax.XTickLabel = xticklabels;
    set(gca,'FontSize',14); 
    plot(x, dfn_U, '-')
    ax.ColorOrderIndex = 1;
    plot(xq, dfn_U_interp, '.--', x, dfe, '-k')% values of numerical difference scheme
    legend('Upwind','Upwind Int','Exact','Location','best')
    title(['Upwind Interpolated n=',num2str(n)]);
    hold off
    
% Downwind 
    subplot(3,2,j+2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 2;
    xlim([0,2*pi])
    set(gca,'XTick',xticks)
    ax.XTickLabel = xticklabels;
    set(gca,'FontSize',14); 
    plot(x, dfn_D, '-')
    ax.ColorOrderIndex = 2;
    plot(xq, dfn_D_interp, '.--', x, dfe, '-k')
    legend('Downwind','Downwind Int','Exact','Location','best')
    title(['Downwind Interpolated n=',num2str(n)]);
    hold off
    
% Central
    subplot(3,2,j+4)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 3;
    xlim([0,2*pi])
    title(['n=',num2str(n)]);
    set(gca,'XTick',xticks)
    ax.XTickLabel = xticklabels;
    set(gca,'FontSize',14); 
    plot(x, dfn_C, '-')
    ax.ColorOrderIndex = 3;
    plot(xq, dfn_C_interp, '.--', x, dfe, '-k')
    legend('Central','Central Int','Exact','Location','best')
    title(['Central Interpolated n=',num2str(n)]);

    
end
    
% Plotting of error over grid spacing in log scale
fig_error = figure;
    loglog(error(:,1),error(:,2),'-',error(:,1),error(:,3),'-',...
        error(:,1),error(:,4),'-',error(:,1),error(:,8),'-')
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    loglog(error(:,1),error(:,5),'--',error(:,1),error(:,6),'--',...
        error(:,1),error(:,7),'--',error(:,1),error(:,9),'--')
%     title('Error plot first derivative');
    xlabel('Grid spacing h');
    ylabel('Relative error');
    legend('Upwind','Downwind','Central','Central Second',...
        'Upwind Int','Downwind Int','Central Int','Central Second Int',...
        'Location','SouthEast')
    
%% Order of error

order_of_er_U = (log10(error(end,2)) - log10(error(end-1,2))) / (log10(error(end,1)) - log10(error(end-1,1)));
order_of_er_D = (log10(error(end,3)) - log10(error(end-1,3))) / (log10(error(end,1)) - log10(error(end-1,1)));
order_of_er_C = (log10(error(end,4)) - log10(error(end-1,4))) / (log10(error(end,1)) - log10(error(end-1,1)));
order_of_er_C2 = (log10(error(end,8)) - log10(error(end-1,8))) / (log10(error(end,1)) - log10(error(end-1,1)));

%% Save files for report

print(fig_first_derivative,'-dpng',"Plots_one/First_derivative.png",'-r150');
print(fig_second_derivative,'-dpng',"Plots_one/Second_derivative.png",'-r150');
print(fig_error,'-dpng',"Plots_one/Error.png",'-r150');
print(fig_first_derivative_int,'-dpng',"Plots_one/fig_first_derivative_int.png",'-r150');



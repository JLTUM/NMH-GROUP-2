function [ run, flow1 ] = time_step_rk1( firstcall, constants, grid, run, ...
        flow1, bconds )
%TIME_STEP_RK Perform one time step with the given explicit Runge-Kutta
%scheme

%% ---- Coefficients of Runge-Kutta scheme --------------------------------
% Butcher tableau of Runge-Kutta scheme by Williamson (1980)
c = [0;1/3;3/4];
a = [0,0,0; 1/3,0,0; -3/16,15/16,0];
b = [1/6,3/10,8/15];

% number of stages
s = length(c);

%% ---- Auxiliary fields --------------------------------------------------
% The time integration method requires some auxiliary 2D fields. To speed
% up the code, we do not want MATLAB to create and destroy these arrays in
% every time step. Therefore, we create persistent variables that remain in
% memory.

% Auxiliary fields for time integration
persistent dhdt;
persistent dhudt;
persistent dhvdt;
persistent hstage;
persistent hustage;
persistent hvstage;

if firstcall
    % Create arrays to store s time derivatives of the different
    % Runge-Kutta substeps
    dhdt  = zeros( [ size(flow1.h),  s] );
    dhudt = zeros( [ size(flow1.hu), s] );
    dhvdt = zeros( [ size(flow1.hv), s] );
    
    % Create arrays to store intermediate values of the 2D fields
    hstage  = zeros( size(flow1.h)  );
    hustage = zeros( size(flow1.hu) );
    hvstage = zeros( size(flow1.hv) );
end

% For safety, we set those fields to NaN in every time step. This should
% prevent bugs where we forget to update the values inside the fields.
dhdt(:,:)  = NaN;
dhudt(:,:) = NaN;
dhvdt(:,:) = NaN;
hstage(:,:) = NaN;
hustage(:,:) = NaN;
hvstage(:,:) = NaN;

%% ---- Loop over substeps ------------------------------------------------
for irk = 1:s
    % compute necessary combinations for current stage
    tstage = run.t + c(irk) * run.dt;
    
    hstage(:,:) = flow1.h;
    hustage(:,:) = flow1.hu;
    hvstage(:,:) = flow1.hv;

    for j = 1:irk-1
        hstage  = hstage + dhdt(:,:,j) * a(irk, j) * run.dt;
        hustage = hustage + dhudt(:,:,j) * a(irk, j) * run.dt;
        hvstage = hvstage + dhvdt(:,:,j) * a(irk, j) * run.dt;
    end

    
    % apply boundary conditions on the intermediate fields
    [hstage, hustage, hvstage] = apply_boundary_conditions(constants, grid, bconds, hstage, hustage, hvstage);
    

    % evaluate stage
    [ dhdt(:,:,irk), dhudt(:,:,irk), dhvdt(:,:,irk)] = ...
        shallow_water_equations( constants, grid, ...
        tstage, hstage, hustage, hvstage, flow1.kst, flow1.zb, ...
        dhdt(:,:,irk), dhudt(:,:,irk), dhvdt(:,:,irk) );
end

%% Combine stages
run.t  = run.t + run.dt;

for j = 1:s
    flow1.h  = flow1.h  + dhdt(:,:,j) * b(j) * run.dt;
    flow1.hu = flow1.hu + dhudt(:,:,j) * b(j) * run.dt;
    flow1.hv = flow1.hv + dhvdt(:,:,j) * b(j) * run.dt;
end

% apply boundary conditions on the final fields
[flow1.h, flow1.hu, flow1.hv] = apply_boundary_conditions(constants, grid, bconds, flow1.h, flow1.hu, flow1.hv);

end

function [ flow1,flow2,flow3 ] = set_initial_condition( grid, flow )
%SET_INITIAL_CONDITION Set initial fields

% ---- Create fields for shallow water equations ----------------------
% Flow depth (cell-centred)
flow1.h = zeros( length(grid.x), length(grid.y) );
flow2.h = zeros( length(grid.x), length(grid.y) );
flow3.h = zeros( length(grid.x), length(grid.y) );
% Specific discharge in x-direction (cell-centred)
flow1.hu = zeros( length(grid.x), length(grid.y) );
flow2.hu = zeros( length(grid.x), length(grid.y) );
flow3.hu = zeros( length(grid.x), length(grid.y) );
% Specific discharge in y-direction (cell-centred)
flow1.hv = zeros( length(grid.x), length(grid.y) );
flow2.hv = zeros( length(grid.x), length(grid.y) );
flow3.hv = zeros( length(grid.x), length(grid.y) );

% Strickler value (cell-centred)
flow1.kst = zeros( length(grid.x) - 2 * grid.NGHOST, length(grid.y) - 2 * grid.NGHOST );
flow2.kst = zeros( length(grid.x) - 2 * grid.NGHOST, length(grid.y) - 2 * grid.NGHOST );
flow3.kst = zeros( length(grid.x) - 2 * grid.NGHOST, length(grid.y) - 2 * grid.NGHOST );
% Bottom elevation (cell-centred)
flow1.zb = zeros( length(grid.x) , length(grid.y) );
flow2.zb = zeros( length(grid.x) , length(grid.y) );
flow3.zb = zeros( length(grid.x) , length(grid.y) );
% Strickler value (cell-centred)
flow1.kst = ones( grid.nx, grid.ny );
flow2.kst = ones( grid.nx, grid.ny );
flow3.kst = ones( grid.nx, grid.ny );
% ---- Fields initialization ----------------------

global infilename

% Case 1
if strcmp(infilename,'infile_2D_swe_damBreak.mat')
 
    disp('Reading Case 1')
 % bottom elevation
%     flow.zb = 
 % kst
    flow.kst(:,:) = 30;
 % water level 
    flow.h(:,:) = 0.1;
    flow.h(1:51,:) = 1;
 % specific discharge
    flow.hu(:,:) = 0;
    
% Case 2
elseif strcmp(infilename,'infile_2D_swe_damBreak_V1.mat')
    
    disp('Reading Case 2-different velocity')
 % bottom elevation
    flow1.I_S = -0.001;
    flow2.I_S = -0.001;
    flow3.I_S = -0.001;
    flow1.zb = repmat(flow1.I_S .* grid.x',1,length(grid.y));
    flow2.zb = repmat(flow2.I_S .* grid.x',1,length(grid.y));
    flow3.zb = repmat(flow3.I_S .* grid.x',1,length(grid.y));
 % kst
    flow1.kst(:,:) = 30;
    flow2.kst(:,:) = 30;
    flow3.kst(:,:) = 30;
 % water level 
    flow1.h(:,:) = 0.01;
    flow1.h(1:51,:) = 1;
    flow2.h(:,:) = 0.01;
    flow2.h(1:51,:) = 0.8;
     flow3.h(:,:) = 0.01;
    flow3.h(1:51,:) = 0.4;
 % specific discharge
    flow1.hu(:,:) = 0;
    flow2.hu(:,:) = 0;
     flow3.hu(:,:) = 0;
   
% Case 3
elseif strcmp(infilename,'infile_2D_swe_damBreak_V2.mat')
    
   disp('Reading Case 3-different kst')
 % bottom elevation
    flow1.I_S = -0.001;
    flow2.I_S = -0.001;
    flow3.I_S = -0.001;
    flow1.zb = repmat(flow1.I_S .* grid.x',1,length(grid.y));
    flow2.zb = repmat(flow2.I_S .* grid.x',1,length(grid.y));
    flow3.zb = repmat(flow3.I_S .* grid.x',1,length(grid.y));
 % kst
    flow1.kst(:,:) = 10;
    flow2.kst(:,:) = 30;
    flow3.kst(:,:) = 110;
 % water level 
    flow1.h(:,:) = 0.1;
    flow1.h(1:51,:) = 1;
    flow2.h(:,:) = 0.1;
    flow2.h(1:51,:) = 1;
    flow3.h(:,:) = 0.1;
    flow3.h(1:51,:) = 1;
   
 % specific discharge
    flow1.hu(:,:) = 0;
    flow2.hu(:,:) = 0;
     flow3.hu(:,:) = 0;
    
% Case 4
elseif strcmp(infilename,'infile_2D_swe_damBreak_V3.mat')
    
    disp('Reading Case 4-different h0')
 % bottom elevation
    flow1.I_S = -0.001;
    flow2.I_S = -0.001;
    flow3.I_S = -0.001;
    flow1.zb = repmat(flow1.I_S .* grid.x',1,length(grid.y));
    flow2.zb = repmat(flow2.I_S .* grid.x',1,length(grid.y));
    flow3.zb = repmat(flow3.I_S .* grid.x',1,length(grid.y));
 % kst
    flow1.kst(:,:) = 30;
    flow2.kst(:,:) = 30;
    flow3.kst(:,:) = 30;
 % water level 
    flow1.h(:,:) = 0.05;
    flow1.h(1:51,:) = 1;
    flow2.h(:,:) = 0.1;
    flow2.h(1:51,:) = 1;
     flow3.h(:,:) = 0.2;
    flow3.h(1:51,:) = 1;
 % specific discharge
    flow1.hu(:,:) = 0;
    flow2.hu(:,:) = 0;
     flow3.hu(:,:) = 0;
   
    

elseif strcmp(infilename,'infile_2D_swe_damBreak_V4.mat')
    
    disp('Reading Case 5')
 % bottom elevation
    flow.I_S = -0.001;
    flow.zb = repmat(flow.I_S .* grid.x',1,length(grid.y));
 % kst
    flow.kst(:,:) = 2;
 % water level 
    flow.h(:,:) = 0.1;
    flow.h(1:51,:) = 1;
 % specific discharge
    flow.hu(:,:) = 0;
   
    elseif strcmp(infilename,'infile_2D_swe_damBreak_V5.mat')
    
    disp('Reading Case 6')
 % bottom elevation
    flow.I_S = -0.001;
    flow.zb = repmat(flow.I_S .* grid.x',1,length(grid.y));
    flow.zb(130:134,:) = 0.09;
 % kst
    flow.kst(:,:) = 30;
 % water level 
    flow.h(:,:) = 0.1;
    flow.h(1:51,:) = 1;
    flow.h(130:134,:) = 0.01;
 
 % specific discharge
    flow.hu(:,:) = 0;
    
        elseif strcmp(infilename,'infile_2D_swe_damBreak_V6.mat')
    
    disp('Reading Case 7')
 % bottom elevation
    flow.I_S = -0.001;
    flow.zb = repmat(flow.I_S .* grid.x',1,length(grid.y));
    flow.zb(130:134,:) = 0.09;
 % kst
    flow.kst(:,:) = 30;
 % water level 
    flow.h(:,:) = 0.1;
    flow.h(1:51,:) = 1;
    flow.h(130:134,:) = 0.01;
 
 % specific discharge
    flow.hu(:,:) = 0;
end

end

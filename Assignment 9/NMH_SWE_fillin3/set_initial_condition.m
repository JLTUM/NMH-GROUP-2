function [ flow ] = set_initial_condition( grid, flow )
%SET_INITIAL_CONDITION Set initial fields

% ---- Create fields for shallow water equations ----------------------
% Flow depth (cell-centred)
flow.h = zeros( length(grid.x), length(grid.y) );

% Specific discharge in x-direction (cell-centred)
flow.hu = zeros( length(grid.x), length(grid.y) );

% Specific discharge in y-direction (cell-centred)
flow.hv = zeros( length(grid.x), length(grid.y) );

% Strickler value (cell-centred)
flow.kst = zeros( length(grid.x) - 2 * grid.NGHOST, length(grid.y) - 2 * grid.NGHOST );

% Bottom elevation (cell-centred)
flow.zb = zeros( length(grid.x) , length(grid.y) );

% Strickler value (cell-centred)
flow.kst = ones( grid.nx, grid.ny );

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
    
    disp('Reading Case 2')
 % bottom elevation
    flow.I_S = -0.001;
    flow.zb = repmat(flow.I_S .* grid.x',1,length(grid.y));
 % kst
    flow.kst(:,:) = 30;
 % water level 
    flow.h(:,:) = 0.1;
    flow.h(1:51,:) = 1;
 % specific discharge
    flow.hu(:,:) = 0;
   
% Case 3
elseif strcmp(infilename,'infile_2D_swe_damBreak_V2.mat')
    
    disp('Reading Case 3')
 % bottom elevation
    flow.I_S = -0.001;
    flow.zb = repmat(flow.I_S .* grid.x',1,length(grid.y));
 % kst
    flow.kst(:,:) = 50;
 % water level 
    flow.h(:,:) = 0.1;
    flow.h(1:51,:) = 1;
 % specific discharge
    flow.hu(:,:) = 0;
    
% Case 4
elseif strcmp(infilename,'infile_2D_swe_damBreak_V3.mat')
    
    disp('Reading Case 4')
 % bottom elevation
    flow.I_S = -0.001;
    flow.zb = repmat(flow.I_S .* grid.x',1,length(grid.y));
 % kst
    flow.kst(:,:) = 80;
 % water level 
    flow.h(:,:) = 0.1;
    flow.h(1:51,:) = 1;
 % specific discharge
    flow.hu(:,:) = 0;
    

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

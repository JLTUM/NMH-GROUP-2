
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
if strcmp(infilename,'infile_2D_swe_channelFlow1.mat')
    flow.zb = repmat(flow.I_S .* grid.x',1,length(grid.y));
    kst = 30;
    flow.kst = kst * flow.kst;
    flow.h(:,:) = 1.35;
    flow.hu(:,:) = 1.5;
    disp('Reading Case 1')
    
% Case 2
elseif strcmp(infilename,'infile_2D_swe_channelFlow2.mat')
    flow.zb = readmatrix('zb2.txt')';
    flow.kst = readmatrix('kst2.txt')';
    flow.I_S = -0.015125;
    flow.h(:,:) = 1.35;
    flow.hu(:,:) = 1.5;
    disp('Reading Case 2')
  
% Case 3
elseif strcmp(infilename,'infile_2D_swe_channelFlow3.mat')
    flow.zb = readmatrix('zb3.txt')';
    flow.kst = readmatrix('kst3.txt')';
    flow.h(:,:) = 0.5;
    flow.hu(:,:) = 1.5;
    disp('Reading Case 3')
    
% Case 4   
elseif strcmp(infilename,'infile_2D_swe_channelFlow4.mat')
    flow.zb = readmatrix('zb4.txt')';
    flow.kst = readmatrix('kst4.txt')';
    flow.h(:,:) = 1.35;
    flow.hu(:,:) = 1.5;
    disp('Reading Case 4')
    
% Case 5
elseif strcmp(infilename,'infile_2D_swe_channelFlow5.mat')
    flow.zb = readmatrix('zb5.txt')';
    flow.kst = readmatrix('kst5.txt')';
    flow.h(:,:) = 1.35;
    flow.hu(:,:) = 1.5;
    disp('Reading Case 5')
    
end


end


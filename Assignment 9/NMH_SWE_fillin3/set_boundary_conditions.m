function [ bconds ] = set_boundary_conditions()
%SET_BOUNDARY_CONDITIONS Set boundary conditions

bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

global infilename

% Case 1
if strcmp(infilename,'infile_2D_swe_damBreak.mat')
    
    bconds.bwest = {'WALL'};
    bconds.beast = {'WALL'};
    
% Case 2
elseif strcmp(infilename,'infile_2D_swe_damBreak_V1.mat')
   
    bconds.bwest = {'WALL'};
    bconds.beast = {'WALL'};
    
%     bconds.bwest = {'HFIX','HUFIX'};
%     bconds.beast = {'WALL'};
%     bconds.hwest = 1.0;
%     bconds.huwest = 0;

% Case 3
elseif strcmp(infilename,'infile_2D_swe_damBreak_V2.mat')
    
    bconds.bwest = {'WALL'};
    bconds.beast = {'WALL'};
    
% Case 4
elseif strcmp(infilename,'infile_2D_swe_damBreak_V3.mat')
    
    bconds.bwest = {'WALL'};
    bconds.beast = {'WALL'};
    
   % Case 5
elseif strcmp(infilename,'infile_2D_swe_damBreak_V4.mat')
    
    bconds.bwest = {'WALL'};
    bconds.beast = {'WALL'};
    
    % Case 6
elseif strcmp(infilename,'infile_2D_swe_damBreak_V5.mat')
    
    bconds.bwest = {'WALL'};
    bconds.beast = {'WALL'};
    
elseif strcmp(infilename,'infile_2D_swe_damBreak_V6.mat')
    
    bconds.bwest = {'WALL'};
    bconds.beast = {'WALL'};
end

end

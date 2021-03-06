function [ bconds ] = set_boundary_conditions()
%SET_BOUNDARY_CONDITIONS Set boundary conditions

bconds.bwest = {'HFIX','HUFIX'};
bconds.beast = {'HEXT','HUEXT'}; 
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

global infilename

% Case 1
if strcmp(infilename,'infile_2D_swe_channelFlow1.mat')
    bconds.hwest = 1.3;
    bconds.huwest = 1.5;
    
% Case 2
elseif strcmp(infilename,'infile_2D_swe_channelFlow2.mat')
    bconds.hwest = 0.351;
    bconds.huwest = 1.5;
    
% Case 3
elseif strcmp(infilename,'infile_2D_swe_channelFlow3.mat')
    bconds.hwest = 1.3;
    bconds.huwest = 1.5;
    
% Case 4
elseif strcmp(infilename,'infile_2D_swe_channelFlow4.mat')
    bconds.hwest = 0.6;
    bconds.huwest = 1.5;
    
% Case 5
elseif strcmp(infilename,'infile_2D_swe_channelFlow5.mat')
    bconds.hwest = 0.4;
    bconds.huwest = 1.5;
end

end


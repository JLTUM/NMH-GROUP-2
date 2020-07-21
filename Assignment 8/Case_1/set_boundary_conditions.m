function [ bconds ] = set_boundary_conditions()
%SET_BOUNDARY_CONDITIONS Set boundary conditions

bconds.bwest = {'HFIX','HUFIX'};
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

global infilename

% Case 1
if strcmp(infilename,'infile_2D_swe_channelFlow1.mat')
    bconds.beast = {'HFIX','HUEXT'}; 
    bconds.hwest = 1.3;
    bconds.huwest = 1.5;
    bconds.heast = 4.1429;
    
% Case 2
elseif strcmp(infilename,'infile_2D_swe_channelFlow2.mat')
    bconds.beast = {'HEXT','HUEXT'}; 
    bconds.hwest = 0.351;
    bconds.huwest = 1.5;
    
% Case 3
elseif strcmp(infilename,'infile_2D_swe_channelFlow3.mat')
    bconds.beast = {'HFIX','HUEXT'}; 
    bconds.hwest = 1.3;
    bconds.huwest = 1.5;
    bconds.heast = 1.3;
    
% Case 4
elseif strcmp(infilename,'infile_2D_swe_channelFlow4.mat')
    bconds.beast = {'HFIX','HUEXT'}; 
    bconds.hwest = 0.6;
    bconds.huwest = 1.5;
    bconds.heast = 0.7;
    
% Case 5
elseif strcmp(infilename,'infile_2D_swe_channelFlow5.mat')
    bconds.beast = {'HEXT','HUEXT'}; 
    bconds.hwest = 0.4;
    bconds.huwest = 1.5;
end

end

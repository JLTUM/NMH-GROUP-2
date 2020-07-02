function [ bconds ] = set_boundary_conditions()
%SET_BOUNDARY_CONDITIONS Set boundary conditions

bconds.bwest = {'HFIX'};
bconds.bwest = {'HUFIX'};    %% HFIX -> h  HUFIX -> hu
bconds.beast = {'HEXT'};
bconds.beast = {'HUEXT'};      %% HEXT   HUEXT extrapolates values
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

 bconds.hwest = 1.3;
 bconds.huwest = 1.5;
% bconds.heast = ???
% bconds.hueast = ???

end


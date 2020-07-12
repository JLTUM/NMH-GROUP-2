function [CFL_max] = compute_CFL_number(constants, grid, dt, h, hu, hv)
%COMPUTE_CFL_NUMBER Computes the CFL number. This is probably perfect.

u = hu ./ h;
v = hv ./ h;
c = sqrt(constants.g .* h);

% compute CFL
CFL = max( abs(u-c), abs(u+c) ) * dt / grid.dx ...
    + max( abs(v-c), abs(v+c) ) * dt / grid.dy;

% use worst case
CFL_max = max(max(CFL));
        
end
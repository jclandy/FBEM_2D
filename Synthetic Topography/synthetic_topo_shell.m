function [x,y,z] = synthetic_topo_shell(op_mode,topo_type,sigma_surf,l_surf,dx)

%% Generates xyz vectors of synthetic sea ice topography

% Virtual sea ice surface topography with Gaussian, lognormal or fractal
% roughness properties

% Input:
% op_mode: 1 = pulse-limited, 2 = SAR
% topo_type: 1 = Gaussian, 2 = lognormal, 3 = fractal
% sigma_surf = rms roughness height
% l_surf = correlation length
% H_surf = Hurst parameter for fractal surfaces
% dx = grid resolution

% Output:
% xyz = surface topography coordinate vertices

% (C) Jack Landy, University of Bristol, 2018

warning('off','all')

%% Create grid

% (Can be modified but these are default grid sizes)
if op_mode == 1
    L = 8000; % across-track diameter of grid, m
    W = 8000; % along-track diameter of grid, m
else
    L = 8000; % across-track diameter of grid, m
    W = 400; % along-track diameter of grid, m
end

[x,y] = meshgrid(-W/2+dx/2:dx:W/2-dx/2,-L/2+dx/2:dx:L/2-dx/2);

%% Generate topography

if topo_type == 1 % Gaussian
    [z,~,~] = rsgene2D_anisotrop(W/dx,L/dx,W,L,sigma_surf,l_surf,l_surf,0);
elseif topo_type == 2 % Lognormal
    [z,~,~] = rsgene2D_anisotrop(W/dx,L/dx,W,L,sigma_surf,l_surf,l_surf,1);
end

warning('on','all')

end


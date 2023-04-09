%SelfCode

%% Create a triangular mesh
% Function createRectangularMesh is used to setup a simple triangular mesh. The
% mesh is visualised in the figure below. Each element (a triangle) and
% boundary element (a line) in the mesh has a unique index that can be
% used to set their properties. The indices of the boundary elements are
% shown in the figure.
%
% <<edge.png>>
%
clc;
clear all;
close all;

xsize = 1;	% width of the region [mm]
ysize = 1;	% height of the region [mm]
dh = 1/2^6;         % discretisation size [mm]
vmcmesh = createRectangularMesh(xsize, ysize, dh);

%% Give optical parameters
% Constant optical parameters are set troughout the medium.

vmcmedium.absorption_coefficient_ex_sol = 0.026;     % absorption coefficient [1/mm]
vmcmedium.absorption_coefficient_ex_f = 0.005;     % absorption coefficient [1/mm]
vmcmedium.absorption_coefficient_em_sol = 0.031;     % absorption coefficient [1/mm]
vmcmedium.scattering_coefficient_ex = 98.4;      % scattering coefficient [1/mm]
vmcmedium.scattering_coefficient_em = 98.4;      % scattering coefficient [1/mm]
vmcmedium.scattering_anisotropy = 0.9;       % anisotropy parameter g of
                                             % the Heneye-Greenstein scattering
                                             % phase function [unitless]
vmcmedium.refractive_index = 1;            % refractive index [unitless]

vmcboundary = createBoundary(vmcmesh,vmcmedium);

%% Create a light source
% Set up a 'cosinic' light source to boundary elements number 4,5,6 and 7.
% This means that the initial propagation direction with respect to the surface
% normal follows a cosine distribution. The photons are launched from random locations
% at these boundary elements.
% line_start = [0 0.65];
% line_end = [0 0.35];
% line_width =1/2^7;
% srcloc = findBoundaries(vmcmesh, 'direction', ...
%                               line_start, ...
%                               line_end,  ...
%                               line_width);

chk1=vmcmesh.r;
chk2=vmcmesh.BH;

lightsource=find(abs(chk1(chk2(:,1),2)-0.5)<=1e-7 & abs(chk1(chk2(:,2),2)-0.5)<=1e-7);
vmcboundary.lightsource(lightsource) = {'direct'};

line_start=[0,0.6];
line_end=[0,0.4];
% Create a direction vector for the light using the line that was used to
% search boundary elements
lightsource_direction = line_end - line_start;

vmcboundary.lightsource_direction(lightsource,1) = lightsource_direction(1);
vmcboundary.lightsource_direction(lightsource,2) = lightsource_direction(2);

% % x-component of the direction
% vmcboundary.lightsource_direction(srcloc,1) = lightsource_direction(1);
% % y-component of the dircetion
% vmcboundary.lightsource_direction(srcloc,2) = lightsource_direction(2);

vmcboundary.lightsource_direction_type(lightsource) = {'absolute'};

%% Set up the boundary
% top_edge_start = [-5,5];
% top_edge_end = [5,5];
% line_width = 0.005;
% top_edge = findBoundaries(vmcmesh, 'direction', top_edge_start, top_edge_end, line_width);
vmcboundary.exterior_refractive_index = 1;

%%
options.frequency=0;    
options.photon_count = 1e6;
options.NBin2Dtheta=16;
% Run the Monte Carlo simulation
% Use the parameters that were generated to run the simulation in the mesh.
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary,options);
%node_solution = nodalbasis2D(vmcmesh,solution);

%% Plot the solution

% rad_full_monte=solution.R_element_fluence;
% flu_full_monte=solution.element_fluence;
% 
% save('Monte_Carlo_comparison_2D_Swapan_version.mat','rad_full_monte','flu_full_monte','vmcmesh')

% Elem2Node = zeros(size(vmcmesh.H,1),size(vmcmesh.r,1));
% for ii=1:size(vmcmesh.H,1)
%     Elem2Node(ii,vmcmesh.H(ii,:)) = 1/3;
% end
% 
% Fluence_nodal = Elem2Node'*solution.element_fluence;
% 
% imagesc(Xvec,fliplr(Yvec),real(reshape(Fluence_nodal,length(Xvec),length(Yvec))));


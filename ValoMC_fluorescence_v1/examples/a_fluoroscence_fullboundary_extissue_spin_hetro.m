%% Working with voxel fomat
% Often imaging data is provided in voxel format.  However, as
% ValoMC uses tetrahedrons as the basis elements, the data are not
% directly compatible. This example demonstrates how to move between
% the two formats.

%% Creating a rectangular 3d mesh
% To create a mesh that can be easily mapped to a voxel grid the
% function createGridMesh can be used

clear all;  
% Create a finer mesh
bins=50.0;
d_size=2.5;
h_disc=d_size/bins;
h_disc_1=d_size/bins;
x_arr = h_disc/2:h_disc:d_size-h_disc/2;
y_arr = h_disc/2:h_disc:d_size-h_disc/2;
z_arr = h_disc_1/2:h_disc_1:d_size-h_disc_1/2;
vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
ck=find(abs(vmcmesh.r(:))<10^-14);
vmcmesh.r(ck)=0.0;

%% Create an anisotropic parameter distribution 
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function

%% Accessing elements using one dimensional indexing
% Note that since there are six times as many tetrahedrons as there are grid
% cells, vmcmedium.absorption_coefficient is six times bigger than F
% A complete assignment can be achieved by repeating the array F six times
mua=0.05;
mus=50.0;
mua_mult=1;
mus_mult=1;
muaf=0.147;
muaf_mult=1;

vmcmedium.scattering_coefficient_ex = mus;
vmcmedium.scattering_coefficient_em = mus*mus_mult;
vmcmedium.absorption_coefficient_ex_sol = mua; % repeat six times
vmcmedium.absorption_coefficient_ex_f = muaf; % repeat six times
%vmcmedium.absorption_coefficient_ex_f =0.0; % repeat six time
vmcmedium.absorption_coefficient_em_sol = mua_mult*mua+muaf*muaf_mult; % repeat six times
%vmcmedium.absorption_coefficient_em_sol = mua_mult*mua; % repeat six times

vmcmedium.scattering_anisotropy = 0.9;        
vmcmedium.refractive_index = 1.33;

vmcboundary = createBoundary(vmcmesh, vmcmedium);   % create a boundary for the mesh
vmcmedium=createMedium(vmcmesh, vmcmedium);
vmcboundary.exterior_refractive_index = 1.33;
% Create a light source
% lightsource = findBoundaries(vmcmesh, 'direction', [0.5 0.5 -0.5], [0.5 0.5 0.5], 1);
% vmcboundary.lightsource(lightsource) = {'direct'};
%% Creating full boundary suorce
z_coord=vmcmesh.r(:,3);
temp=vmcmesh.BH;
z_coord_1=z_coord(temp(:,1));
z_coord_2=z_coord(temp(:,2));
z_coord_3=z_coord(temp(:,3));

chk_1=find((z_coord_1==0) & (z_coord_2==0) & (z_coord_3==0));
lightsource=chk_1;

% chk_2=z_coord(temp(chk_1,:));
vmcboundary.lightsource(lightsource) = {'direct'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
options.frequency=0;    
options.photon_count = 6e8;
options.Qyield_f = 1;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary,options);
% node_solution=nodalbasis3D(vmcmesh,solution);
% nod_solution=nodalbasis(vmcmesh,solution);

%% ploting new way
nvoxels= length(vmcmesh.H) / 6;
nx=int16(nvoxels^(1/3));
ny=int16(nvoxels^(1/3));
nz=int16(nvoxels^(1/3));
first = reshape(solution.element_fluence(1:nvoxels), nx,ny,nz);
second = reshape(solution.element_fluence(nvoxels+1:2*nvoxels),nx,ny,nz);
third = reshape(solution.element_fluence(nvoxels*2+1:3*nvoxels),nx,ny,nz);
fourth = reshape(solution.element_fluence(nvoxels*3+1:4*nvoxels),nx,ny,nz);
fifth = reshape(solution.element_fluence(nvoxels*4+1:5*nvoxels),nx,ny,nz);
sixth = reshape(solution.element_fluence(nvoxels*5+1:6*nvoxels),nx,ny,nz);
solution.grid_fluence = (first+second+third+fourth+fifth+sixth)/6;
%solution.grid_fluence = solution.grid_fluence1/max(max(max(solution.grid_fluence1)));
fn=reshape(solution.grid_fluence(int16(bins/2),:,:),size(solution.grid_fluence,1),size(solution.grid_fluence,3));
figure(3);clf
imagesc(z_arr',x_arr', log10(fn));
%J_im1 = imrotate(J_im,90,'bilinear','crop');
xlabel('z [mm]');
ylabel('x [mm]');
%figure
%imshow(J_im1)
colorbar
%colormap(makec2f)

% %slice(X, Y, Z, solution.grid_fluence,  [0.5156], [], []);
% xlabel('x [mm]');
% ylabel('y [mm]');
% zlabel('z [mm]');
% 
% view(125,25);
% snapnow;
%% for flouroscence
first1 = reshape(solution.F_element_fluence(1:nvoxels), nx,ny,nz);
second2 = reshape(solution.F_element_fluence(nvoxels+1:2*nvoxels),nx,ny,nz);
third3 = reshape(solution.F_element_fluence(nvoxels*2+1:3*nvoxels),nx,ny,nz);
fourth4 = reshape(solution.F_element_fluence(nvoxels*3+1:4*nvoxels),nx,ny,nz);
fifth5 = reshape(solution.F_element_fluence(nvoxels*4+1:5*nvoxels),nx,ny,nz);
sixth6 = reshape(solution.F_element_fluence(nvoxels*5+1:6*nvoxels),nx,ny,nz);
solution.grid_fluence1 = (first1+second2+third3+fourth4+fifth5+sixth6)/6;
%solution.grid_fluence1 = solution.grid_fluence1/max(max(max(solution.grid_fluence1)));
%fn2=reshape(sum(solution.grid_fluence1,2),size(solution.grid_fluence1,1),size(solution.grid_fluence1,3));
fn1=reshape(solution.grid_fluence1(int16(bins/2),:,:),size(solution.grid_fluence1,1),size(solution.grid_fluence1,3));
figure(4);clf
imagesc(z_arr',x_arr', log10(fn1));
%J_im1 = imrotate(J_im,90,'bilinear','crop');
xlabel('z [mm]');
ylabel('x [mm]');
%figure
%imshow(J_im1)
colorbar
%clim([-0.5 2])
%colormap(makec2f)

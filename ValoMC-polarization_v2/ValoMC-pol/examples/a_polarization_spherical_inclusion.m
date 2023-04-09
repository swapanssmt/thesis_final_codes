%% Working with voxel fomat
% Often imaging data is provided in voxel format.  However, as
% ValoMC uses tetrahedrons as the basis elements, the data are not
% directly compatible. This example demonstrates how to move between
% the two formats.

%% Creating a rectangular 3d mesh
% To create a mesh that can be easily mapped to a voxel grid the
% function createGridMesh can be used

clear all;  
%% test
bins=100.0;
h_disc=1.0/bins;
h_disc_1=0.12/bins;
x_arr = h_disc/2:h_disc:1.0-h_disc/2;
y_arr = h_disc/2:h_disc:1.0-h_disc/2;
z_arr = h_disc_1/2:h_disc_1:0.12-h_disc_1/2;
vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
ck=find(abs(vmcmesh.r(:))<10^-14);
vmcmesh.r(ck)=0.0;
%% Calculating parameter by Mie Theory
% background
radius=double(0.05);
lambda=double(0.6328);
nre_med=double(1.33);
nim_med=double(0.0);
nre_p=double(1.59);
nim_p=double(0.0);
nangles=int64(1000);
rho=double(19.11);
[mus_1, s11_1, s12_1, s33_1, s43_1] = MTmex(radius, lambda, nre_med, nim_med, nre_p, nim_p, nangles, rho);
% spherical inclusion
radius_2=double(1.0);
rho_2=double(1.11e-3);
[mus_2, s11_2, s12_2, s33_2, s43_2] = MTmex(radius_2, lambda, nre_med, nim_med, nre_p, nim_p, nangles, rho_2);
% s11_2=zeros(nangles+1,1);
% s12_2=zeros(nangles+1,1);
% s33_2=zeros(nangles+1,1);
% s43_2=zeros(nangles+1,1);
s11=real([s11_1 s11_2]);
s12=real([s12_1 s12_2]);
s33=real([s33_1 s33_2]);
s43=real([s43_1 s43_2]);
%% Create an anisotropic parameter distribution 
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function

%% Accessing elements using one dimensional indexing
% Note that since there are six times as many tetrahedrons as there are grid
% cells, vmcmedium.absorption_coefficient is six times bigger than F
% A complete assignment can be achieved by repeating the array F six times
vmcmedium.scattering_coefficient = mus_1;
vmcmedium.absorption_coefficient = 0.05; % repeat six times
vmcmedium.scattering_anisotropy = 0.9;
vmcmedium.layer = 1;
vmcmedium.refractive_index = 1.33;

%%
vmcboundary = createBoundary(vmcmesh, vmcmedium);   % create a boundary for the mesh
vmcmedium=createMedium(vmcmesh, vmcmedium);
vmcboundary.exterior_refractive_index = 1.33;
%% creating spherical inclusion
% radius_s=0.05;
% center_s=[0.6 0.6 0.06];
% sphere=findElements(vmcmesh, 'sphere', center_s,radius_s);
% vmcmedium.scattering_coefficient(sphere) = mus_2;
% vmcmedium.layer(sphere) = 2;
%%
% Create a light source
% lightsource = findBoundaries(vmcmesh, 'direction', [0.5 0.5 -0.5], [0.5 0.5 0.5], 1);
% vmcboundary.lightsource(lightsource) = {'direct'};

z_coord=vmcmesh.r(:,3);
temp=vmcmesh.BH;
z_coord_1=z_coord(temp(:,1));
z_coord_2=z_coord(temp(:,2));
z_coord_3=z_coord(temp(:,3));

ck_1=find((z_coord_1 == 0) & (z_coord_2 == 0) & (z_coord_3 == 0));
lightsource=ck_1;
% chk_2=z_coord(temp(chk_1,:));
vmcboundary.lightsource(lightsource) = {'direct'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Creating full boundary suorce
% z_coord=vmcmesh.r(:,3);
% temp=vmcmesh.BH;
% z_coord_1=z_coord(temp(:,1));
% z_coord_2=z_coord(temp(:,2));
% z_coord_3=z_coord(temp(:,3));
% 
% chk_1=find((z_coord_1==0) & (z_coord_2==0) & (z_coord_3==0));
% lightsource=chk_1;
% 
% % chk_2=z_coord(temp(chk_1,:));
% vmcboundary.lightsource(lightsource) = {'direct'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
options.s11=s11;
options.s12=s12;
options.s33=s33;
options.s43=s43;
options.s0=[1.0; 1.0; 0.0; 0.0];
options.nangles=nangles;
options.activate_pol=1; % use 0 for running without polarization and 1 for running with polarization
options.frequency=0;    
options.photon_count = 2e8;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary,options);
%% z-boundary
z_coord=vmcmesh.r(:,3);
temp=vmcmesh.BH;
z_coord_1=z_coord(temp(:,1));
z_coord_2=z_coord(temp(:,2));
z_coord_3=z_coord(temp(:,3));

chk_1=find((z_coord_1==0) & (z_coord_2==0) & (z_coord_3==0));
z_pboundary=chk_1;
%% calculating sum
RIB=sum(solution.IB(chk_1))/options.photon_count;
RQB=sum(solution.QB(chk_1))/options.photon_count;
RUB=sum(solution.UB(chk_1))/options.photon_count;
RVB=sum(solution.VB(chk_1))/options.photon_count;
chk_2=find((z_coord_1==5) & (z_coord_2==5) & (z_coord_3==5));
TIB=sum(solution.IB(chk_2))/options.photon_count;
TQB=sum(solution.QB(chk_2))/options.photon_count;
TUB=sum(solution.UB(chk_2))/options.photon_count;
TVB=sum(solution.VB(chk_2))/options.photon_count;
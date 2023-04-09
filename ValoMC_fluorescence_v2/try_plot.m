% cla
% long = vmcmesh.r(:,1)';                              % longitude data
% lat =  vmcmesh.r(:,2)';                               % latitude data
% rural = vmcmesh.r(:,3)';                     % percent rural data
% fatalities =node_solution.element_fluence;% fatalities data
% 
% scatter3(long,lat,rural,60,fatalities,'filled')    % draw the scatter plot
% 
% ax = gca;
% ax.XDir = 'reverse';
% %view(-31,14)
% xlabel('W. Longitude')
% ylabel('N. Latitude')
% zlabel('% Rural Population')
% 
% cb = colorbar;                                     % create and label the colorbar
% cb.Label.String = 'Fatalities per 100M vehicle-miles';

% xslice = [0.2 0.5 0.7 0.9];                               % define the cross sections to view
% yslice = [];
% zslice = ([1 0]);
% % x = vmcmesh.r(:,1)';                              % longitude data
% % y =  vmcmesh.r(:,2)';                               % latitude data
% % z = vmcmesh.r(:,3)';                     % percent rural data
% temp =node_solution.element_fluence; 
% 
% slice(vmcmesh.r, temp, xslice, yslice, zslice)    % display the slices
% %ylim([-3 3])
% view(-34,24)
% 
% cb = colorbar;                                  % create and label the colorbar
% cb.Label.String = 'Temperature, C';

% x = -0.5:;
% y =vmcmesh.r(:,2)';
% z = vmcmesh.r(:,3)';
% [x,y,z] = meshgrid(x,y,z) ; 
% f = x.*exp(x.^2 + y.^2 + z.^2);
% figure
% hold on
% for i = 1:size(x,3)
%     surf(x(:,:,i),y(:,:,i),z(:,:,i),f(:,:,i))
% end
% view(3)
% shading interp

%%
h_disc=1/2^4;
x_arr = h_disc/2:h_disc:1-h_disc/2;
y_arr = h_disc/2:h_disc:1-h_disc/2;
z_arr = h_disc/2:h_disc:1-h_disc/2;
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function
grid_fluence = repmat(node_solution.element_fluence,,size(X));
slice(X, Y, Z, grid_fluence, [0 0.2 0.5 0.7 1], [], 1);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');

view(125,25);
snapnow;
function elements = findsphere(vmcmesh,center,radius)
nvoxels= length(vmcmesh.H) / 6;
elements1=[];
for el=1:nvoxels
 first =  el;
 second = el + nvoxels;
 third =  el + nvoxels*2;
 fourth = el + nvoxels*3;
 fifth =  el + nvoxels*4;
 sixth =  el + nvoxels*5;
 voxel1=[vmcmesh.H(first,:) vmcmesh.H(second,:) vmcmesh.H(third,:) vmcmesh.H(fourth,:) vmcmesh.H(fifth,:) vmcmesh.H(sixth,:)]';
 voxel=unique(voxel1);
 if(size(voxel,1)~=8)
     fprintf('something wrong \n');
 end
 x_cood=(vmcmesh.r(voxel(1),1)+vmcmesh.r(voxel(2),1)+vmcmesh.r(voxel(3),1)+vmcmesh.r(voxel(4),1)+vmcmesh.r(voxel(5),1)+vmcmesh.r(voxel(6),1)+vmcmesh.r(voxel(7),1)+vmcmesh.r(voxel(8),1))/8;
 y_cood=(vmcmesh.r(voxel(1),2)+vmcmesh.r(voxel(2),2)+vmcmesh.r(voxel(3),2)+vmcmesh.r(voxel(4),2)+vmcmesh.r(voxel(5),2)+vmcmesh.r(voxel(6),2)+vmcmesh.r(voxel(7),2)+vmcmesh.r(voxel(8),2))/8;
 z_cood=(vmcmesh.r(voxel(1),3)+vmcmesh.r(voxel(2),3)+vmcmesh.r(voxel(3),3)+vmcmesh.r(voxel(4),3)+vmcmesh.r(voxel(5),3)+vmcmesh.r(voxel(6),3)+vmcmesh.r(voxel(7),3)+vmcmesh.r(voxel(8),3))/8;
 centroid_voxel=[x_cood y_cood z_cood];
 dist1=centroid_voxel-center;
 dist=sqrt(sum(dist1.^2));
 
 if (dist<=radius)
     elements1(end+1)=first;
     elements1(end+1)=second;
     elements1(end+1)=third;
     elements1(end+1)=fourth;
     elements1(end+1)=fifth;
     elements1(end+1)=sixth;

 end
end
elements=elements1';
end
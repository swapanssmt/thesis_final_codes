radius_s=2;
center_s=[12.75 12.75 18.75];
a_chk=findElements(vmcmesh, 'sphere', center_s,radius_s);
plot3([vmcmesh.r(vmcmesh.H(a_chk,1),1) vmcmesh.r(vmcmesh.H(a_chk,2),1) vmcmesh.r(vmcmesh.H(a_chk,3),1)  vmcmesh.r(vmcmesh.H(a_chk,4),1)]', ...
    [vmcmesh.r(vmcmesh.H(a_chk,1),2) vmcmesh.r(vmcmesh.H(a_chk,2),2) vmcmesh.r(vmcmesh.H(a_chk,3),2) vmcmesh.r(vmcmesh.H(a_chk,4),2)]', ...
    [vmcmesh.r(vmcmesh.H(a_chk,1),3) vmcmesh.r(vmcmesh.H(a_chk,2),3) vmcmesh.r(vmcmesh.H(a_chk,3),3) vmcmesh.r(vmcmesh.H(a_chk,4),3)]','o')
hold on
plot3(center_s(1),center_s(2),center_s(3),'o-')
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
hold off
function node_solution = nodalbasis2D(vmcmesh,solution)
sz=size(vmcmesh.r,1);
asize=size(solution.R_element_fluence,2);
node_solution.element_fluence=zeros(sz,1);
node_solution.R_element_fluence=zeros(sz,asize);
node_solution.boundary_exitance=zeros(sz,1);
node_solution.R_boundary_exitance=zeros(sz,asize);
node_solution.boundary_fluence=zeros(sz,1);
node_solution.R_boundary_fluence=zeros(sz,asize);
%size(unique(vmcmesh.BH),1);
% size(node_solution.element_fluence)
% size(node_solution.R_element_fluence)
for i=1:size(solution.element_fluence,1)
    a1=solution.element_fluence(i,:)/3;
    node_solution.element_fluence(vmcmesh.H(i,1),:)=node_solution.element_fluence(vmcmesh.H(i,1),:)+a1;
    node_solution.element_fluence(vmcmesh.H(i,2),:)=node_solution.element_fluence(vmcmesh.H(i,2),:)+a1;
    node_solution.element_fluence(vmcmesh.H(i,3),:)=node_solution.element_fluence(vmcmesh.H(i,3),:)+a1;
    a2=solution.R_element_fluence(i,:)/3;
    node_solution.R_element_fluence(vmcmesh.H(i,1),:)=node_solution.R_element_fluence(vmcmesh.H(i,1),:)+a2;
    node_solution.R_element_fluence(vmcmesh.H(i,2),:)=node_solution.R_element_fluence(vmcmesh.H(i,2),:)+a2;
    node_solution.R_element_fluence(vmcmesh.H(i,3),:)=node_solution.R_element_fluence(vmcmesh.H(i,3),:)+a2;
end
for i=1:size(solution.boundary_exitance,1)
    a3=solution.boundary_exitance(i,:)/2;
    node_solution.boundary_exitance(vmcmesh.BH(i,1),:)=node_solution.boundary_exitance(vmcmesh.BH(i,1),:)+a3;
    node_solution.boundary_exitance(vmcmesh.BH(i,2),:)=node_solution.boundary_exitance(vmcmesh.BH(i,2),:)+a3;
    a4=solution.R_boundary_exitance(i,:)/2;
    node_solution.R_boundary_exitance(vmcmesh.BH(i,1),:)=node_solution.R_boundary_exitance(vmcmesh.BH(i,1),:)+a4;
    node_solution.R_boundary_exitance(vmcmesh.BH(i,2),:)=node_solution.R_boundary_exitance(vmcmesh.BH(i,2),:)+a4;
    a5=solution.boundary_fluence(i,:)/2;
    node_solution.boundary_fluence(vmcmesh.BH(i,1),:)=node_solution.boundary_fluence(vmcmesh.BH(i,1),:)+a5;
    node_solution.boundary_fluence(vmcmesh.BH(i,2),:)=node_solution.boundary_fluence(vmcmesh.BH(i,2),:)+a5;
    a6=solution.R_boundary_fluence(i,:)/2;
    node_solution.R_boundary_fluence(vmcmesh.BH(i,1),:)=node_solution.R_boundary_fluence(vmcmesh.BH(i,1),:)+a6;
    node_solution.R_boundary_fluence(vmcmesh.BH(i,2),:)=node_solution.R_boundary_fluence(vmcmesh.BH(i,2),:)+a6;
end
end
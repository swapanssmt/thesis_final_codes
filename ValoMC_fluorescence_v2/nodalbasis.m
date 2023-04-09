function nod_solution = nodalbasis(vmcmesh,solution)
sz=size(vmcmesh.r,1);
asize=size(solution.element_radiance,2);
nod_solution.element_fluence=zeros(sz,1);
nod_solution.element_radiance=zeros(sz,asize);
nod_solution.boundary_exitance=zeros(sz,1);
nod_solution.R_boundary_exitance=zeros(sz,asize);
nod_solution.boundary_fluence=zeros(sz,1);
nod_solution.boundary_radiance=zeros(sz,asize);
nod_solution.F_element_fluence=zeros(sz,1);
nod_solution.F_element_radiance=zeros(sz,asize);
nod_solution.F_boundary_exitance=zeros(sz,1);
nod_solution.F_R_boundary_exitance=zeros(sz,asize);
nod_solution.F_boundary_fluence=zeros(sz,1);
nod_solution.F_boundary_radiance=zeros(sz,asize);
for i=1:sz
    elements=[];
    for j=1:size(vmcmesh.H,2)
        elements=union(elements,find(vmcmesh.H(:,j)==i));
    end
    nod_solution.element_fluence(i)=sum(solution.element_fluence(elements))/size(elements,1);
    nod_solution.element_radiance(i,:)=sum(solution.element_radiance(elements))/size(elements,1);
    nod_solution.F_element_fluence(i)=sum(solution.F_element_fluence(elements))/size(elements,1);
    nod_solution.F_element_radiance(i,:)=sum(solution.F_element_radiance(elements))/size(elements,1);
    elements_b=[];
    for k=1:size(vmcmesh.BH,2)
        elements_b=union(elements_b,find(vmcmesh.BH(:,k)==i));
    end
    nod_solution.boundary_exitance(i)=sum(solution.boundary_exitance(elements_b))/size(elements_b,1);
    nod_solution.R_boundary_exitance(i,:)=sum(solution.R_boundary_exitance(elements_b))/size(elements_b,1);
    nod_solution.boundary_fluence(i)=sum(solution.boundary_fluence(elements_b))/size(elements_b,1);
    nod_solution.boundary_radiance(i,:)=sum(solution.boundary_radiance(elements_b))/size(elements_b,1);
    nod_solution.F_boundary_exitance(i)=sum(solution.F_boundary_exitance(elements_b))/size(elements_b,1);
    nod_solution.F_R_boundary_exitance(i,:)=sum(solution.F_R_boundary_exitance(elements_b))/size(elements_b,1);
    nod_solution.F_boundary_fluence(i)=sum(solution.F_boundary_fluence(elements_b))/size(elements_b,1);
    nod_solution.F_boundary_radiance(i,:)=sum(solution.F_boundary_radiance(elements_b))/size(elements_b,1);
end
end

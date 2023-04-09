function centroid = angle3Ddiscreteization()
load('discreteization_to_compare.mat');
pa;
ta;
x1=pa(ta(:,1),1);
x2=pa(ta(:,2),1);
x3=pa(ta(:,3),1);
y1=pa(ta(:,1),2);
y2=pa(ta(:,2),2);
y3=pa(ta(:,3),2);
z1=pa(ta(:,1),3);
z2=pa(ta(:,2),3);
z3=pa(ta(:,3),3);
%centroid1=(x1+x2+x3)/3;
centroid1=[(x1+x2+x3)/3,(y1+y2+y3)/3,(z1+z2+z3)/3];
centroid=centroid1;
for i=1:length(centroid)
    nm=norm(centroid(i,:));
    centroid(i,:)=centroid(i,:)/nm;
end
end
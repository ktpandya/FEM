%% Mesh Creation
clear all
clc
close all

model=createpde;
importGeometry(model,'squaredomain.stl');
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.05);

E1=findElements(mesh,'region','Face',1);
N1=findNodes(mesh,'region','Face',1);
figure
pdemesh(mesh,'ElementLabels','on','NodeLabels','on')
hold on
pdemesh(mesh.Nodes(:,N1),mesh.Elements(:,E1))
[p,e,t]=meshToPet(mesh);
coordinates=p;
connectivity=t(1:3,:);
nNodes=size(coordinates,2);
nElements=size(t,2);
a=0;

for i=1:nNodes
 if(p(1,i)==0||p(1,i)==1)
    a=a+1;
    nodesonboundary(1,a)=i;


   
 else
    if(p(2,i)==0||p(2,i)==1)
        a=a+1;
        nodesonboundary(1,a)=i;
        
    end
 end
end
disp(a)
nBoundaryNodes=a;
fileID=fopen('coordinates.txt','w');
fprintf(fileID,'%12.8f %12.8f\n',p);
fclose(fileID);
fileID=fopen('connectivity.txt','w');
fprintf(fileID,'%d %d %d\n',t(1:3,:));
fclose(fileID);
fileID=fopen('boundarynodes.txt','w');
fprintf(fileID,'%d\n',nodesonboundary);
fclose(fileID);


















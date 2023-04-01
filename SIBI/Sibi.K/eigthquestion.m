% clear all;
% clc;
% close all;
% syms x1 y1 x2 y2 x3 y3 x y X Y;
% f(X,Y)=X+Y-X^2-Y^2;
% phix=(x1-x3)*x + (x2-x3)*y + x3;
% phiy=(y1-y3)*x+(y2-y3)*y +y3;
% h=subs(f,X,phix);
% fnew=subs(h,Y,phiy);
% N0= x;
% N1= y;
% N2= 1-x-y;
% h0=fnew*N0;
% h1=fnew*N1;
% h2=fnew*N2;
% f1=int(int(h0,x,0,1-y),y,0,1);
% f2=int(int(h1,x,0,1-y),y,0,1);
% f3=int(int(h2,x,0,1-y),y,0,1);
% disp("done");

model=createpde;
importGeometry(model,'squaredomain.stl');
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.05);

E1=findElements(mesh,'region','Face',1);
N1=findNodes(mesh,'region','Face',1);
figure
pdemesh(mesh,'ElementLabels','on','NodeLabels','on')
hold on
pdemesh(mesh.Nodes(:,N1),mesh.Elements(:,E1))

uh=importdata("sol.dat");
figure(1)
pdeplot(model,"XYData",uh)
figure(2)
pdeplot(model,"ZData",uh)

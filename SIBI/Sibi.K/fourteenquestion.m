clear all;
clc;
close all;
Hmax=zeros(1,10);
RMSERROR=zeros(1,10);
Hmax(1,1)=1;
RMSERROR(1,1)=0.000931680083505064;
Hmax(1,2)=0.5;
RMSERROR(1,2)=0.000399301052415516;
Hmax(1,3)=0.2;
RMSERROR(1,3)=0.000204150821802324;
Hmax(1,4)=0.1;
RMSERROR(1,4)=4.13438495401547e-05;
Hmax(1,5)=0.05;
RMSERROR(1,5)=8.71793313056647e-06;

loglog(Hmax,RMSERROR)
set ( gca, 'xdir', 'reverse' )
xlabel("log(Hmax) (Decreasing →)")
ylabel("log(RMSERROR)")
title("log(RMSERROR) Vs log(Hmax)")
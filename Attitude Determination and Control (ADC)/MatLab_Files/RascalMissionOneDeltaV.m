%RascalMissionOneDeltaV
%Rev -
%Last Edited: 3/3/2014, Tom Moline
%Functions Called: CWPrussingRendezvousSolver.m

%Desription================================================================
%The CWPrussingRendezvvousSolverTester script assigns varying levels of
%initial and final displacement between two spacecraft and plots the result
%of varying the initial relative velocity between tehm
%==========================================================================

%Begin Code

%Clear Screen and Variables
clear
clc

%Initialize Variables
rtgt_=[6697.4756; 1794.5831; 0.0];
drFinalMin=0.01;
drFinalMax=1;
dv0Max=0.01;
Nr=5;
Nv=20;
Nt=200;
tmax=100;

dr0_=linspace(drFinalMin,drFinalMax,Nr);
dr0Matrix=[dr0_;dr0_;dr0_];

figure (1)
CWSolverESC(drFinalMin,drFinalMax,dv0Max,Nr,Nv,Nt,rtgt_,tmax)
figure (2)
CWSolverSK(dr0Matrix,Nt,rtgt_,tmax)
figure (3)
CWSolverRDZ(drFinalMax,drFinalMin,Nr,Nt,rtgt_,tmax)

% t=linspace(.1*60,tmax*60,Nt);
% dr0_=linspace(drFinalMin,drFinalMax,Nr);
% dr0Matrix=[dr0_;dr0_;dr0_];
% dr0mag=sqrt(sum(abs(dr0Matrix).^2));
% dv0_=[0;0;0];
% dvmag=zeros(length(dr0_),length(t));
% drmag=zeros(length(dr0_),length(t));
% 
% for i=1:length(dr0_)
%     for j=1:length(t)
%     [dr_, dv_]=CWPrussing(dr0Matrix(:,i),dv0_,rtgt_,t(j));
%     drmag(i,j)=(dr_(1)^2+dr_(2)^2+dr_(3)^2)^.5;
%     dvmag(i,j)=(dv_(1)^2+dv_(2)^2+dv_(3)^2)^.5;
%     end
% end
% 
% figure (2)
% surf(t./60,dr0mag,drmag,'FaceAlpha',0.6)









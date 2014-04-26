function [deltaVSKError]=CWSolverSKError(dr0,dv0,Nt,rtgt_,tmax)
%The CWSolverSK function takes in information on the maximum initial
%separtion relative velocity and desired relative final positoin of two s/c
%and, based on the time over which one wants to simulate, plots Nr cases of
%total deltaV plots to station keep at varous final relative positons.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      drFinalMin Min Final Rel Pos            Scalar             km
%      drFinalMax Max Final Rel Pos            Scalar             km
%      dv0Max     Max Initial Rel Vel          Scalar             km/s
%      Nr         Number of Position Cases     Scalar           Unitless
%      Nv         Number of Velocity Cases     Scalar           Unitless
%      Nt         Number of Time Cases         Scalar           Unitless
%      rtgt_      Inertial Target Position     3x1 Vector         km
%      tmax       Simulation Time              Sclarl             mins

%Initial Release, CWSolverSK.m, Tom Moline, 3/3/2014

%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
dv0mag=sqrt(sum(abs(dv0).^2));
drf=dr0;
t=linspace(60*.1,60*tmax,Nt); %Simultation time, s
deltaVSKError=zeros(length(dv0mag),length(t)); %Pre-Allocate for speed

for i=1:length(dv0mag)
    for j=1:length(t)
        %Call CWPrussingRendezvousSolver Function
        [deltaVi,deltaVf]=CWPrussingRendezvousSolver(dr0,dv0(:,i),rtgt_,drf,t(j));
        deltaVSKError(i,j)=sum(sqrt((abs(deltaVi)+abs(deltaVf)).^2)); %Find deltaV for each case
        if deltaVSKError(i,j)>=.005 %Ignore values greater than 0.150 km/s
            deltaVSKError(i,j)=.005;
        end
    end
end

function CWSolverRDZ(drInitialMin,drInitialMax,Nr,Nt,rtgt_,tmax)
%The CWSolverSK function takes in information on the maximum initial
%separtion relative velocity and desired relative final positoin of two s/c
%and, based on the time over which one wants to simulate, plots Nr cases of
%total deltaV plots to station keep at varous final relative positons.
%
%==========================================================================
% Variable Name     Variable Description      Variable Type  Variable Units
%==========================================================================
%      drInitialMin Min Initial Rel Pos          Scalar             km
%      drFinalMax   Max Initial Rel Pos          Scalar             km
%      dv0Max       Max Initial Rel Vel          Scalar             km/s
%      Nr           Number of Position Cases     Scalar           Unitless
%      Nv           Number of Velocity Cases     Scalar           Unitless
%      Nt           Number of Time Cases         Scalar           Unitless
%      rtgt_        Inertial Target Position     3x1 Vector         km
%      tmax         Simulation Time              Sclarl             mins

%Initial Release, CWSolverRDZ.m, Tom Moline, 3/3/2014

%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
dr0=linspace(drInitialMin,drInitialMax,Nr); %Final separtion vector, km
dr0matrix=[dr0;dr0;dr0]; %Matrix of final separation vectors
dr0mag=sqrt(sum(abs(dr0matrix).^2));
dv0_=[0;0;0];
drf_=[0;0;0]; %Case for s/c starting out together
t=linspace(60*.1,60*tmax,Nt); %Simultation time, s
deltaVRDZ=zeros(length(dr0mag),length(t)); %Pre-Allocate for speed

%Fill deltV structure field with diffrent values for each drfx case within
%the drfmatrix.
for k=1:length(dr0mag)
        for j=1:length(t)
            %Call CWPrussingRendezvousSolver Function
            [deltaViRDZ,deltaVfRDZ]=CWPrussingRendezvousSolver(dr0matrix(:,k),dv0_,rtgt_,drf_,t(j));
            deltaVRDZ(k,j)=sum(sqrt((abs(deltaViRDZ)+abs(deltaVfRDZ)).^2)); %Find deltaV for each case
        if deltaVRDZ(k,j)>=.05 %Ignore values greater than 0.150 km/s
            deltaVRDZ(k,j)=.05;
        end
        end
end %Move to next drfx value

%Plot Results
surf(t./60,dr0mag*1000,deltaVRDZ*1000,'FaceAlpha',0.6)
grid on
set(gca,'GridLineStyle','-')
xlabel('Transfer Time (min)')
ylabel('Initial Relative Displacement Magnitude (m)')
zlabel('Total Delta V (m/s)')
s=sprintf('Total DeltaV Required for Rendezvous with Initial Relative Displacements Between %.2f m and %.2f m',1000*dr0mag(1),...
    1000*dr0mag(length(dr0mag)));
title(s)
colorbar
view(18,26)

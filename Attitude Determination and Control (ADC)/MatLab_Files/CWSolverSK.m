function CWSolverSK(dr0,Nt,rtgt_,tmax)
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
dv0_=[0;0;0];
dr0mag=sqrt(sum(abs(dr0).^2));
drf=dr0;
t=linspace(60*.1,60*tmax,Nt); %Simultation time, s
deltaVSK=zeros(length(dr0mag),length(dr0mag)); %Pre-Allocate for speed

for i=1:length(dr0mag)
    for j=1:length(t)
        %Call CWPrussingRendezvousSolver Function
        [deltaVi,deltaVf]=CWPrussingRendezvousSolver(dr0(:,i),dv0_,rtgt_,drf(:,i),t(j));
        deltaVSK(i,j)=sum(sqrt((abs(deltaVi)+abs(deltaVf)).^2)); %Find deltaV for each case
        if deltaVSK(i,j)>=.01 %Ignore values greater than 0.150 km/s
            deltaVSK(i,j)=.01;
        end
    end
end

%Plot Results
surf(t./60,dr0mag*1000,deltaVSK*1000,'FaceAlpha',0.6)
grid on
set(gca,'GridLineStyle','-')
xlabel('Transfer Time (min)')
ylabel('Initial Relative Position Magnitude (m)')
zlabel('Total Delta V (m/s)')
s=sprintf('Total DeltaV Required for Station Keeping at Initial Relative Displacements between %.2f m and %.2f m',1000*dr0mag(1),...
    1000*dr0mag(length(dr0mag)));
title(s)
view(18,26)


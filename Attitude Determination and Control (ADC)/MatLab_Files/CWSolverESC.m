function CWSolverESC(drFinalMin,drFinalMax,dv0Max,Nr,Nv,Nt,rtgt_,tmax)
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
drfx=linspace(drFinalMin,drFinalMax,Nr); %Final separtion vector, km
drfmatrix=[drfx;drfx;drfx]; %Matrix of final separation vectors
field1='deltaVtotESC'; %Structre field for total deltaV cases
dvx0=linspace(0.001,dv0Max,Nv); %Initial reltative velocity vector, km/s
dvmatrix0=[dvx0;dvx0;dvx0]; %Initial reltative velocity matrix, km/s
dvmag=sqrt(sum(abs(dvmatrix0).^2));
dr0_=[0;0;0]; %Case for s/c starting out together
t=linspace(60*.1,60*tmax,Nt); %Simultation time, s
deltaV=struct(field1,{}); %Create total deltaV structure
deltaVvoxESC=zeros(length(dvx0),length(drfx)); %Pre-Allocate for speed

%Fill deltV structure field with diffrent values for each drfx case within
%the drfmatrix.
for k=1:length(drfx)
    for i=1:length(dvx0)
        for j=1:length(t)
            %Call CWPrussingRendezvousSolver Function
            [deltaViESC,deltaVfESC]=CWPrussingRendezvousSolver(dr0_,dvmatrix0(:,i),rtgt_,drfmatrix(:,k),t(j));
            deltaVvoxESC(i,j)=sum(sqrt((abs(deltaViESC)+abs(deltaVfESC)).^2)); %Find deltaV for each case
        if deltaVvoxESC(i,j)>=.15 %Ignore values greater than 0.150 km/s
            deltaVvoxESC(i,j)=.15;
        end
        end
    end
    deltaV(k).deltaVtotESC=deltaVvoxESC(:,:); %Assign deltaV to strcuture
end %Move to next drfx value

%Plot Results
for i=1:length(drfx)
    hold on;
    surf(t./60,dvmag*1000,deltaV(i).deltaVtotESC*1000,'FaceAlpha',0.6)
    legendinfo{i}=['Final Relative Position (m): ' num2str(drfx(i)*1000)];
end
hold off
grid on
set(gca,'GridLineStyle','-')
xlabel('Transfer Time (min)')
ylabel('Initial Relative Velocity Magnitude (m/s)')
zlabel('Total Delta V (m/s)')
s=sprintf('Total DeltaV Required for Separaton Relative Velociity Magnitudes Between %.2f m/s and %.2f m/s',1000*dvmag(1),...
    1000*dvmag(length(dvx0)));
title(s)
colorbar
legend(legendinfo)
view(18,26)

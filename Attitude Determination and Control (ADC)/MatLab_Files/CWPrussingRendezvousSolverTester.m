%CWPrussingRendezvousSolverTester
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
drfx=linspace(.01,10,5); %Final separtion vector, km
drfmatrix=[drfx;drfx;drfx]; %Matrix of final separation vectors
field1='deltaVtotSK'; %Structre field for total deltaV cases
dvx0=linspace(0.001,0.01,20); %Initial reltative velocity vector, km/s
dvmatrix0=[dvx0;dvx0;dvx0]; %Initial reltative velocity matrix, km/s
dvmag=sqrt(sum(abs(dvmatrix0).^2));
dr0_=[0;0;0]; %Case for s/c starting out together
rtgt_=[6697.4756; 1794.5831; 0.0]; %Target s/c positon, km
t=linspace(60*.1,60*300,200); %Simultation time, s
deltaV=struct(field1,{}); %Create total deltaV structure
deltaVvoxSK=zeros(length(dvx0),length(drfx)); %Pre-Allocate for speed

%Fill deltV structure field with diffrent values for each drfx case within
%the drfmatrix.
for k=1:length(drfx)
    for i=1:length(dvx0)
        for j=1:length(t)
            %Call CWPrussingRendezvousSolver Function
            [deltaViSK,deltaVfSK]=CWPrussingRendezvousSolver(dr0_,dvmatrix0(:,i),rtgt_,drfmatrix(:,k),t(j));
            deltaVvoxSK(i,j)=sum(sqrt((abs(deltaViSK)+abs(deltaVfSK)).^2)); %Find deltaV for each case
        if deltaVvoxSK(i,j)>=.15 %Ignore values greater than 0.150 km/s
            deltaVvoxSK(i,j)=.15;
        end
        end
    end
    deltaV(k).deltaVtotSK=deltaVvoxSK(:,:); %Assign deltaV to strcuture
end %Move to next drfx value

%Plot Results
for i=1:length(drfx)
    hold on;
    surf(t./60,dvmag*1000,deltaV(i).deltaVtotSK*1000,'FaceAlpha',0.6)
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






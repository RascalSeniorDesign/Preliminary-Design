%% Script Description
% The IdealRascalMission script produces the deltaV's, relative postions,
% and realtive velocities associated with executing the ideal Rascal
% mission.

%% Functions Used
% CWSolverIdeal.m
% CWPrussingSolver.m
% CWPrussingDeltaVSolver.m
% CWStateMatrix.m

%% Revision History
% Revision Number: -
% Revision Comment: Initial Release
% Date: 3/16/2014
% Last Edited by: Tom Moline

%% Initialize Variables
clear
clc
clf
% Initial Separation Relative Velocity (km/s)
dv0_=[0.0005;0.0005;0];
% Manevuer Time (mins)
t=90;
% Number of Data Points
Nt=200;
% Target Spacecraft Inertial Position (km)
rtgt0_=[6697.4756; 1794.5831; 0.0];
% Desired Inspection Stationkeeping Distance (km)
drf_=[0.001;0.01;0.001];
% Iniital Relative Position (km)
dr0_=[0;0;0];

%% Create Mission Phase Structure
field1='deltaVi';
field2='deltaVf';
field3='dr_';
field4='dv_';
IdealMission=struct(field1,{},field2,{},field3,{},field4,{});

%% Fill Structure Fields
for i=1:5
    if i==2
        drf_=[0.001;0.01;0.001]; %ISK
    elseif i==3
        drf_=[0.001;.1;0.01]; %Continued Separation
    elseif i==5
        drf_=[0.001;0.01;0.001]; %Rendezvous
    end
  [IdealMission(i).deltaVi,...
 IdealMission(i).deltaVf,...
 IdealMission(i).dr_,...
 IdealMission(i).dv_        ]=CWSolverIdeal(drf_,dv0_,dr0_,t,Nt,rtgt0_); 
  dv0_=[0;0;0]; %Set initial realtive velocity to zero
  dr0_=IdealMission(i).dr_(:,Nt); %Reset initial position
end

%% Print Total DeltaV Used
for i=1:5
    VfViDiff=abs(IdealMission(i).deltaVf-IdealMission(i).deltaVi);
    totaldeltaVm(i)=sqrt(sum(VfViDiff.^2));
end
totaldeltaV=sum(totaldeltaVm)*1000;
fprintf('The Total DeltaV Used for this mission is: %0.2f m/s\n',totaldeltaV)
%% Plot Results
figure(1)
hold on
cc=jet(5);
 for i=1:5
    h=plot3(IdealMission(i).dr_(1,:)*1000,...
         IdealMission(i).dr_(2,:)*1000,...
         IdealMission(i).dr_(3,:)*1000,...
         'color',        cc(i,:),...
         'LineWidth',    2                  );
     k=plot3(IdealMission(i).dr_(1,1)*1000,...
         IdealMission(i).dr_(2,1)*1000,...
         IdealMission(i).dr_(3,1)*1000,...
         'ok',...
         'MarkerSize',      5,...
         'MarkerFaceColor', 'k');
 end
 s=sprintf('Rascal Mission Orbit Path, with Total Delta-V of %0.2f m/s',totaldeltaV);
 xlabel('Relative Cross-Track Displacement (m)')
 ylabel('Relative In-Track Displacement (m)')
 zlabel('Relative Out-of-Plane Dispacement (m)')
 title(s)
 grid on
 set(gca,'GridLineStyle','-')
 view(-40,20)
 legend(k,'Transition Points')
 
 %% Combine Separate Relative Positions into Single Matrix
 dr_=horzcat(IdealMission(:).dr_);
 
 



    

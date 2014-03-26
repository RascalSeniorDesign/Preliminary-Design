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
dv0_=[0;0.001;0];
% Manevuer Time (mins)
t=90;
% Number of Data Points
Nt=200;
% Target Spacecraft Inertial Position (km)
rtgt0_=[6697.4756; 1794.5831; 0.0];
% Iniital Relative Position (km)
dr0_=[0;0;0];
% Mission Type (5=No Docking, 6=Docking)
Nmission=6;
% Pertubation (km)
Pert=0.005;
% Desired Inspection Stationkeeping Distance (km)
drf_=[Pert;0.01;Pert];

%% Create Mission Phase Structure
field1='deltaVi';
field2='deltaVf';
field3='dr_';
field4='dv_';
IdealMission=struct(field1,{},field2,{},field3,{},field4,{});

%% Fill Structure Fields
for i=1:Nmission
    if i==2
        drf_=[Pert;0.01;Pert]; %ISK
        dv0_=[0;0.0001;0];
    elseif i==3
        drf_=[Pert;.1;Pert]; %Continued Separation
    elseif i==5
        drf_=[Pert;0.01;Pert]; %Rendezvous
    elseif i>5
        drf_=[0;0;0];
    end
  [IdealMission(i).deltaVi,...
 IdealMission(i).deltaVf,...
 IdealMission(i).dr_,...
 IdealMission(i).dv_        ]=CWSolverIdeal(drf_,dv0_,dr0_,t,Nt,rtgt0_); 
  dv0_=[0;0;0]; %Set initial realtive velocity to zero
  dr0_=IdealMission(i).dr_(:,Nt); %Reset initial position
end

%% Print Total DeltaV Used for Each Case
deltaVi=horzcat(IdealMission(:).deltaVi);
deltaVf=horzcat(IdealMission(:).deltaVf);
deltaVcomp=abs(deltaVi)+abs(deltaVf);
deltaVmag=sqrt(sum(deltaVcomp.^2))*1000;
fprintf('The tolta DeltaV used during Initial Separation is: %0.4f m/s\n',deltaVmag(1))
fprintf('The tolta DeltaV used during ISK is: %0.4f m/s\n',deltaVmag(2))
fprintf('The tolta DeltaV used during Continued Separation is: %0.4f m/s\n',deltaVmag(3))
fprintf('The tolta DeltaV used during RSK is: %0.4f m/s\n',deltaVmag(4))
fprintf('The tolta DeltaV used during Rendezvous is: %0.4f m/s\n',deltaVmag(5))
if Nmission==6
    fprintf('The tolta DeltaV used during Docking is: %0.4f m/s\n\n',deltaVmag(6))
end

if Nmission==5
    deltaVtotalmission=sum(deltaVmag)+deltaVmag(2)*28+2*sum(deltaVmag(2:5));
else
    deltaVtotalmission=2*sum(deltaVmag)+14*deltaVmag(2)+sum(deltaVmag(2:5));
end
fprintf('Then the total DeltaV for the entire mission is: %0.4f m/s\n\n',deltaVtotalmission)
%% Plot Results
figure(1)
hold on
cc=jet(Nmission);
 for i=1:Nmission
    h=plot3(IdealMission(i).dr_(1,:)*1000,...
         IdealMission(i).dr_(2,:)*1000,...
         IdealMission(i).dr_(3,:)*1000,...
         'color',        cc(i,:),...
         'LineWidth',    3                  );
     k=plot3(IdealMission(i).dr_(1,1)*1000,...
         IdealMission(i).dr_(2,1)*1000,...
         IdealMission(i).dr_(3,1)*1000,...
         'ok',...
         'MarkerSize',      5,...
         'MarkerFaceColor', 'k');
 end
 if Nmission==5
    s=sprintf('Rascal Mission Orbit Path for Phase 1 of CONOPS-1, with Total Delta-V of %0.2f m/s',deltaVtotalmission);
 else
    s=sprintf('Rascal Mission Orbit Path for Phase 1 of CONOPS-2, with Total Delta-V of %0.2f m/s',deltaVtotalmission);
 end
 xlabel('Relative Cross-Track Displacement (m)')
 ylabel('Relative In-Track Displacement (m)')
 zlabel('Relative Out-of-Plane Dispacement (m)')
 title(s)
 grid on
 set(gca,'GridLineStyle','-')
 view(-40,20)
 

 
 %% Combine Separate Relative Positions into Single Matrix
 dr_=horzcat(IdealMission(:).dr_)*1000;
 dv_=horzcat(IdealMission(:).dv_)*1000;
 t=linspace(0.1*60,(90*Nmission)*60,Nmission*Nt);
 A=vertcat(t,dr_,dv_);
 file=fopen('rascal_plot.txt','w');
 fprintf(file,'stk.v.8.0\r\n\r\nBEGIN Ephemeris\r\n\r\n');
 fprintf(file,'NumberofEphemerisPoints\t%s\r\n',num2str(Nt*Nmission));
 fprintf(file,'ScenarioEpoch\t\t\t20 Oct 2014 13:21:15.0000\r\n');
 fprintf(file,'InterpolationMethod\t\tLagrange\r\n');
 fprintf(file,'InterpolationOrder\t\t5\r\n');
 fprintf(file,'DistanceUnit\t\t\t\tMeters\r\n');
 fprintf(file,'CentralBody\t\t\t\tEarth\r\n');
 fprintf(file,'CoordinateSystem\t\tCW\r\n');
 fprintf(file,'CoordianteSystemEpoch\t20 Oct 2014 13:21:15.0000\r\n\r\n');
 fprintf(file,'EphemerisTimePosVel\r\n\r\n');
 fprintf(file,'%-0.6f\t\t%0.6f\t\t%0.6f\t\t%0.6f\t\t%0.6f\t\t%0.6f\t\t%0.6f\r\n',A);
 fprintf(file,'\r\nEND Ephemeris');
 
 
 



    

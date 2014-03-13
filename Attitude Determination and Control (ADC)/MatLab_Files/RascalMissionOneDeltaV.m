%RascalMissionOneDeltaV
%Rev 1
%Last Edited: 3/12/2014, Tom Moline
%Edit Comments
%1: Added ability to find min and max deltaV for each maneuver case
%Functions Called: CWSolverESC.m, CWSolverRDZ.m, CWSolverSK.m, CWPrussing.m

%Desription================================================================
%The RascalMissionOneDeltaV script calculates the DeltaV associated with
%each of the phases of the Rascal mission, as defined in RCL-O-CMQA3 Rascal
%CONOPS Trade Study. It accomplishes this through the use of linear orbit
%analysis and for the ranges laid out in the variable initialization
%section.
%==========================================================================

%Begin Code

%Clear Screen and Variables
clear
clc

%% Initialize Variables
rtgt_=[6697.4756; 1794.5831; 0.0]; %Target Orbit, km
drFinalMin=0.01; %First Final Displacement Value, km
drFinalMax=.1;%Last Final Displacement Value, km
dv0Max=0.002; %Maximum Iniital Relative Velocity Consideration, km/s
Nr=5; %Number of Displacement Cases
Nv=5; %Number of Velocity Cases
Nt=200; %Number of Transfer Time Cases
tmax=100; %Maximum Transfer Time, mins
t=linspace(0.1,tmax,Nt);

dr0_=linspace(drFinalMin,drFinalMax,Nr); %Initial Relative Displacement
dr0Matrix=[dr0_;dr0_;dr0_]; %Matrix Version

%% Find DeltaV Values Associeated with Escape
figure (1)
deltaV=CWSolverESC(drFinalMin,drFinalMax,dv0Max,Nr,Nv,Nt,rtgt_,tmax);

%Find Reasonable Times at Which to Perform Maneuvers
deltaVESC1=deltaV(1).deltaVtotESC(2,:).*(deltaV(1).deltaVtotESC(Nv,:)<0.0035);
deltaVESC2=deltaV(Nr).deltaVtotESC(1,:).*(deltaV(Nr).deltaVtotESC(1,:)<0.0035);
[deltaVESC1min,tESC1min]=min(deltaV(1).deltaVtotESC(2,(0.8*Nt:Nt)));
[deltaVESC1max,tESC1max]=max(deltaVESC1);
[deltaVESC2min,tESC2min]=min(deltaV(Nr).deltaVtotESC(2,(0.8*Nt:Nt)));
[deltaVESC2max,tESC2max]=max(deltaVESC2);
%% Find DeltaV Values Associated with Stationkeeping
figure (2)
deltaVSK=CWSolverSK(dr0Matrix,Nt,rtgt_,tmax);

%Find Reasonable Times at Which to Perform Maneuvers
deltaVSKtot=deltaVSK.*(deltaVSK<=0.0005);
[deltaVISKmin,tISKmin]=min(deltaVSK(1,(0.8*Nt:Nt)));
[deltaVISKmax,tISKmax]=max(deltaVSKtot(1,:));
[deltaVRSKmin,tRSKmin]=min(deltaVSK(Nr,(0.8*Nt:Nt)));
[deltaVRSKmax,tRSKmax]=max(deltaVSKtot(Nr,:));
%% Find DeltaV Values Associated with Rendezvous
figure (3)
deltaVtot=CWSolverRDZ(drFinalMin,drFinalMax,Nr,Nt,rtgt_,tmax);

%Find Reasonable Times at Which to Perform Maneuvers
deltaVRDZtot=deltaVtot.*(deltaVtot<=0.0025);
[deltaVRDZmax,tRDZmax]=max(deltaVRDZtot(Nr,:));
[deltaVRDZmin,tRDZmin]=min(deltaVtot(Nr,(0.8*Nt:Nt)));

%% Define Manevuer DeltaV Ranges
disp('Each Maneuver for The Rascal Mission Involves the following deltaVs')
fprintf('Initial Separtion to %0.1f m: \nMin: %0.5f m/s for transfer time of %0.1f minutes \nMax: %0.5f m/s for transfer time of %0.1f minutes \n\n',...
    drFinalMin*1000,deltaVESC1min*1000,t(tESC1min+.8*Nt),deltaVESC1max*1000,t(tESC1max))
fprintf('Continued Separtion from %0.1f m to %0.1f m: \nMin: %0.5f m/s for transfer time of %0.1f minutes \nMax: %0.5f m/s for transfer time of %0.1f minutes \n\n',...
    drFinalMin*1000,drFinalMax*1000,deltaVESC2min*1000,t(tESC2min+.8*Nt),deltaVESC2max*1000,t(tESC2max))
fprintf('Inspection Stationkeeping at %0.1f m: \nMin: %0.5f m/s for transfer time of %0.1f minutes \nMax: %0.5f m/s for transfer time of %0.1f minutes \n\n',...
    drFinalMin*1000,deltaVISKmin*1000,t(tISKmin+.8*Nt),deltaVISKmax*1000,t(tISKmax))
fprintf('Remote Stationkeeping at %0.1f m: \nMin: %0.5f m/s for transfer time of %0.1f minutes \nMax: %0.5f m/s for transfer time of %0.1f minutes \n\n',...
    drFinalMax*1000,deltaVRSKmin*1000,t(tRSKmin+.8*Nt),deltaVRSKmax*1000,t(tRSKmax))
fprintf('Rendezvous from %0.1f m: \nMin: %0.5f m/s for transfer time of %0.1f minutes \nMax: %0.5f m/s for transfer time of %0.1f minutes \n\n',...
    drFinalMax*1000,deltaVRDZmin*1000,t(tRDZmin+.8*Nt),deltaVRDZmax*1000,t(tRDZmax))











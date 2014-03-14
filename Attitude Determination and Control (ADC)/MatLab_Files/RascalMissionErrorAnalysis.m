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
clf
%% Initialize Variables

%Define Target Orbit and Initial Conditions
rtgt_=[6697.4756; 1794.5831; 0.0]; %km
dr0=[0.01; 0.01; 0.01]; %km
dv0min=0; %Minimum Initial Relative Velocity, km/s
dv0max=0.001; %Maximum Initial Relative Velocity km/s
Nv=200; %Number of Cases between Initial and Final Value
dv0range=linspace(dv0min,dv0max,Nv);
null=zeros(1,Nv);
dv0=[dv0range;dv0range;dv0range]; %Intial Relative Velocity Matrix
dv0mag=sqrt(sum(abs(dv0).^2)); %Velocit Magnitude Array
tmax=100; %Max Transfer time, mins
Nt=200; %Number of Trnasfer time cases
t=linspace(.1,tmax,Nt);
field1='deltaV';
SKError=struct(field1,{}); %Structure of TotalDeltaV's

%% Fill SKError Structure
for i=1:Nv
    [SKError(i).deltaV]=CWSolverSKError(dr0,dv0(:,i),Nt,rtgt_,tmax);
end

%% Plot Results

%Create Colormap
cc=jet(Nv);
%Create Plot
hold on
figure (1)
for i=1:Nv
    plot(t,SKError(i).deltaV*1000,'color',cc(i,:),'LineWidth',2)
end
%Label the Plot
xlabel('Transfer Time (mins)')
ylabel('Total \DeltaV (m/s)')
%Create Title String that Automatically Updates
s=sprintf(['DeltaV Required for ISK with',...
    'Initial Relative Velocity Magnitudes Between',...
    ' %0.2f m/s and %0.2f m/s \n&\n An In-Track Displacement of',...
    ' %0.1f m, Cross-Track Displacement of %0.1f m, and',...
    ' Out of Plane Displacement of %0.1f m'],...
    dv0mag(1)*1000,dv0mag(Nv)*1000,dr0(2)*1000,dr0(1)*1000,dr0(3)*1000);

set(gca, 'CLim', [dv0mag(1)*1000, dv0mag(Nv)]*1000)
colormap(cc); %Redefine colormap in terms of one being used
c=colorbar; %Add a colorbar
d=ylabel(c,'Initial Relative Velocity Magnitudes (m/s)','rot',-90);
e=get(d,'position'); %Get ylabel position data for colorbar
e(1,1)=e(1,1)+2; %Shift ylabel x-position for colorbar
set(d,'position',e) %Reassign ylabel positon
title(s) %Add title
%% Find Min and Max Value Ranges
minRange=(0.8*Nt:Nt);
maxRange=(0.6*Nt:Nt);
[ISKErrorMinBest, tISKErrorMinBest]=min(SKError(1).deltaV(minRange)*1000);
[ISKErrorMaxBest, tISKErrorMaxBest]=max(SKError(1).deltaV(maxRange)*1000);
[ISKErrorMinWorst, tISKErrorMinWorst]=min(SKError(Nv).deltaV(minRange)*1000);
[ISKErrorMaxWorst, tISKErrorMaxWorst]=max(SKError(Nv).deltaV(maxRange)*1000);

%% Plot and Print the Results

%Plot of Minimum DeltaV Values
plot([t(tISKErrorMinBest+.8*Nt),t(tISKErrorMinWorst+.8*Nt)],...
    [ISKErrorMinBest,ISKErrorMinWorst],'k--o','LineWidth',1,...
    'MarkerFaceColor','k')

%Annotate the Minimum DeltaV Line
text(t(tISKErrorMinWorst+0.8*Nt),ISKErrorMinWorst+.02,...
    'Min DeltaV Line','VerticalAlignment','bottom',...
    'HorizontalAlignment','Right','FontSize',10)

%Print Min and Max DeltaV Values to the Command Window
fprintf(['The DeltaV Values Associated with',...
    'Ideal Initial Conditions are: \n\n'])

fprintf(['Min: %0.2f m/s for transfer time of %0.2f mins\nMax:',...
    '%0.2f m/s for transfer time of %0.2f mins\n\n'],...
    ISKErrorMinBest,t(tISKErrorMinBest+.8*Nt),...
    ISKErrorMaxBest,t(tISKErrorMaxBest+.6*Nt))

fprintf(['The DeltaV Values Associated with',...
    'Worst Initial Conditions are: \n\n'])

fprintf(['Min: %0.2f m/s for transfer time of',...
    '%0.2f mins\nMax: %0.2f m/s for transfer time of %0.2f mins\n\n'],...
    ISKErrorMinWorst,t(tISKErrorMinWorst+.8*Nt),...
    ISKErrorMaxWorst,t(tISKErrorMaxWorst+.6*Nt))







%% RascalMissionErrorAnalysis Script Information
%Last Edited: 3/13/2014, Tom Moline
%Rev: -
%Functions Called: CWSolverSKError.m

%Desription================================================================
%The RascalMissionErrorAnalysis script varies the initial relative
%velocities between a target and interceptor spacecraft at the beginning
%of either an Inspection Stationkeeping (ISK) or Remote Stationkeeping
%(RSK) maneuver, as to assess the sensitivity of said maneuvers to s/c
%positioning error. The output is a plot of detlaV's, as well as a print
%out of the max and min deltaV's associated with best and worst case
%initial position and velocity scenarios.
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
dv0=[dv0range;dv0range;null]; %Intial Relative Velocity Matrix
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
%% Find Min and Max Value Ranges
minRange=(0.8*Nt:Nt);
maxRange=(0.6*Nt:Nt);
[ISKErrorMinBest, tISKErrorMinBest]=min(SKError(1).deltaV(minRange)*1000);
[ISKErrorMaxBest, tISKErrorMaxBest]=max(SKError(1).deltaV(maxRange)*1000);
[ISKErrorMinWorst, tISKErrorMinWorst]=min(SKError(Nv).deltaV(minRange)*1000);
[ISKErrorMaxWorst, tISKErrorMaxWorst]=max(SKError(Nv).deltaV(maxRange)*1000);

%% Print Max and Min Values

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
%% Plot Results
%Create Colormap
cc=jet(Nv);
%Create Plot
hold on
figure (1)
for i=1:Nv
    plot(t,SKError(i).deltaV*1000,'color',cc(i,:),'LineWidth',2)
end
%Plot of Minimum DeltaV Values
plot([t(tISKErrorMinBest+.8*Nt),t(tISKErrorMinWorst+.8*Nt)],...
    [ISKErrorMinBest,ISKErrorMinWorst],'k--o','LineWidth',1,...
    'MarkerFaceColor','k')

hold off
%% Label the Plot

%Crate Min DeltaV Line Annotation
text(t(tISKErrorMinWorst+0.8*Nt),ISKErrorMinWorst+.02,...
    'Min DeltaV Line','VerticalAlignment','bottom',...
    'HorizontalAlignment','Right','FontSize',10)

xlabel('Transfer Time (mins)')
ylabel('Total \DeltaV (m/s)')
%Create Title String that Automatically Updates
s=sprintf(['DeltaV Required for ISK with',...
    'Initial Relative Velocity Magnitudes Between',...
    ' %0.2f m/s and %0.2f m/s \n&\n An In-Track Displacement of',...
    ' %0.1f m, Cross-Track Displacement of %0.1f m, and',...
    ' Out of Plane Displacement of %0.1f m'],...
    dv0mag(1)*1000,dv0mag(Nv)*1000,dr0(2)*1000,dr0(1)*1000,dr0(3)*1000);
title(s) %Add title

%Create and Label the ColorBar
colormap(cc); %Redefine colormap in terms of one being used
c=colorbar; %Add a colorbar
set(gca, 'CLim', [dv0mag(1)*1000, dv0mag(Nv)]*1000) %Set colorbar limits
d=ylabel(c,'Initial Relative Velocity Magnitudes (m/s)','rot',-90);
e=get(d,'position'); %Get ylabel position data for colorbar
e(1,1)=e(1,1)+2; %Shift ylabel x-position for colorbar
set(d,'position',e) %Reassign ylabel positon








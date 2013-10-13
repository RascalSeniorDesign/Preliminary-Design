%=======================Rascal Senior Design===============================
%============================RCL-C-NAV1====================================
%=======================Author: Tom Moline=================================
%=====================Approval: Nate Richard===============================
%=====================Revision: 1 =========================================
%===============Last Edit Date: October 13, 2013===========================

%============================Summary=======================================
%This script serves as a preliminary method of calculating the delta V
%associated with orbital maneuvers made relative to another satellite,
%specifically, that between Butch (a released satellite) and Sundance (a
%secured satellite).

%============================Inputs========================================
% Rp = Orbital Periapsis Altitude (km)
% Ra = Orbital Apoapsis Altitude (km)
% rstar_ =Inital Location of Chase Satellite Relative to Focal Point(km)
% r_ = Initial Location of Target Satellite Relatvie to Focal Point (km)
%    NOTE: This is assumed to be a 3-D Vector input [x, y]
%
% All variables with star in name are of chase satellite. All others are of
% the target satellite.
%
% All variables with an underscore are vectors


%===========================Outputs========================================
% Total Delta V to go from intial positon to one within 10 m of Sundance

%==========================Begin Code======================================

clear
clc

% Create Delta V Structure

deltaV_total=struct('rstar_',{},'dr_',{},'deltav_x',{},...
    'deltav_y',{},'Rp',{},'Ra',{},'period',{},'n',{},'t',{});

Ra=linspace(200,900,10);
Rp=linspace(200,900,10);
dr_=[linspace(.001,1,10);linspace(.001,1,10);linspace(0,0,10)];
rstar_=[linspace(200,900,10);linspace(0,0,10);linspace(0,0,10)];

% Rp = input('Please input the periapsis of a selected orbit: ');
% Ra = input('Please input the apoapsis of a selected orbit: ');

% Calculate Orbital Parameters

mu = 398600;

for i=1:length(Ra)
    e=0;
    a = Ra(i)/(1+e); %Semi-Major Axis, km
    period_value=2*pi*sqrt(a^3/mu);
    deltaV_total(i).period=period_value;
    deltaV_total(i).dr_=dr_(:,i);
    deltaV_total(i).rstar_=rstar_(:,i);
    deltaV_total(i).n=sqrt(mu/(rstar_(1,i))^3);
    deltaV_total(i).Ra=Ra(i);
    deltaV_total(i).Rp=Rp(i);
end


%Indicate desired orbital time step for final rendezvous

tvalue=linspace(0.1,12*pi,250);


for j=1:length(Ra)
    for i=1:length(tvalue)
	
	%Calculate Delta V
    
        deltaV_total(i).t=tvalue(i);

        s=sin(deltaV_total(j).n*tvalue(i)*(180/pi));
        c=cos(deltaV_total(j).n*tvalue(i)*(180/pi));

        M=[(4-3*c) 0 0;6*(s-deltaV_total(j).n*tvalue(i)) 1 0;0 0 c];
        N=[s/deltaV_total(j).n (2/deltaV_total(j).n)*(1-c) 0; -(2/deltaV_total(j).n)*(1-c) (4*s-3*deltaV_total(j).n*tvalue(i))/deltaV_total(j).n 0; 0 0 s/deltaV_total(j).n];
        S=[3*deltaV_total(j).n*s 0 0; -6*deltaV_total(j).n*(1-c) 0 0; 0 0 -deltaV_total(j).n*s];
        T=[c 2*s 0; -2*s 4*c-3 0; 0 0 c];

        deltav=(T/N*M-S)*deltaV_total(:,j).rstar_;
        deltavx(j,i)=deltav(1);
        deltavy(j,i)=deltav(2);
    end
end

for j=1:length(Ra)
    for i=1:length(tvalue)
        deltaV_total(j).deltav_x=deltavx(j,i);
        deltaV_total(j).deltav_y=deltavy(j,i);
    end
end

% for i=1:length(tvalue)
%     plot(deltaV_total(i).t,deltaV_total(i).deltav_x,deltaV_total(i).t,deltaV_total(i).deltav_y)
% end
% xlabel('Time (Period Step Size)')
% ylabel('Delta V (km/s)')


















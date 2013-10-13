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
    'deltav_y',{},'Rp',{},'Ra',{},'period',{},'n',{},'t',{},'deltav',{});

norbits=3;

Ra=linspace(200,900,norbits);
Rp=linspace(200,900,norbits);
dr_=[linspace(.1,.1,norbits);linspace(.1,.1,norbits);linspace(0,0,norbits)];
rstar_=[linspace(200,900,norbits);linspace(0,0,norbits);linspace(0,0,norbits)];

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

        s=sin(deltaV_total(j).n*tvalue(i)*(180/pi));
        c=cos(deltaV_total(j).n*tvalue(i)*(180/pi));

        M=[(4-3*c) 0 0;6*(s-deltaV_total(j).n*tvalue(i)) 1 0;0 0 c];
        N=[s/deltaV_total(j).n (2/deltaV_total(j).n)*(1-c) 0; -(2/deltaV_total(j).n)*(1-c) (4*s-3*deltaV_total(j).n*tvalue(i))/deltaV_total(j).n 0; 0 0 s/deltaV_total(j).n];
        S=[3*deltaV_total(j).n*s 0 0; -6*deltaV_total(j).n*(1-c) 0 0; 0 0 -deltaV_total(j).n*s];
        T=[c 2*s 0; -2*s 4*c-3 0; 0 0 c];

        deltavtemp=(T/N*M-S)*deltaV_total(:,j).rstar_;
        if deltavtemp(1)>1000
            deltavx(j,i)=1000;
        elseif deltavtemp(1)<-1000
            deltavx(j,i)=-1000;    
        else
            deltavx(j,i)=deltavtemp(1);
        end
        
        if deltavtemp(2)>1000
             deltavy(j,i)=1000;
        elseif deltavtemp(2)<-1000
             deltavy(j,i)=-1000;    
        else
             deltavy(j,i)=deltavtemp(2);
        end
    end
    deltaV_total(j).t=tvalue(1,:);
    deltaV_total(j).deltav_x=deltavx(j,:);
    deltaV_total(j).dletav_y=deltavy(j,:);
end

cc=jet(20);

hold on
for i=1:norbits
    plot(tvalue,deltaV_total(i).deltav_x,'color',cc(i,:))
    legendinfo{i}=['Orbit (km): ' int2str(rstar_(1,i))];
end
legend(legendinfo)
xlabel('Time (Period Step Size)')
ylabel('Delta V in Local X Direction (km/s)')

    


















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
% dr_ =Inital displacement vector between each satellite (km)
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
    'deltav_y',{},'Rp',{},'Ra',{},'period',{},'n',{},'t',{},'deltav',{},'numorbits',{},'deltav_total',{},'deltav_z',{},...
    'deltav_xi',{},'deltav_xf',{},'deltav_yi',{},'deltav_yf',{},'deltav_zi',{},'deltav_zf',{});

%Define number of cases to observe (norbits) and Number of Orbits to Plot
%(orbitnumber)

orbitnumber=1;
norbits=1;

%Create initial displacement and velocity vector arrays for each case

r_ = [3212.59    4572    -3877.23]; % km
v_ = [-6.379    1.003    -4.106]; % km/s
dr_=[0 0 0]; %km
dv_=[0 .00050 0];

%Call orbits_plot function to obtain relative distance and velocity data
%between chaser and target satellite for orbit duration

[dr_realtime,dv_realtime,x,y,z,dx,dy,dz,dx_realtime,dy_realtime,dz_realtime,x_realtime,y_realtime,z_realtime,theta,T,e] = orbits_plot(r_,v_,dr_,dv_,orbitnumber);

%Define initial conditions for initial displacement and velocity

dr_=[linspace(mean(x_realtime),mean(x_realtime),norbits);linspace(mean(y_realtime),mean(y_realtime),norbits);linspace(mean(z_realtime),mean(z_realtime),norbits)];
rstar_=[linspace(3212.59,3212.59,norbits);linspace(4572,4572,norbits);linspace(-3877.23,-3877.23,norbits)];
dv_initial=[linspace(mean(dx_realtime),mean(dx_realtime),norbits);linspace(mean(dy_realtime),mean(dy_realtime),norbits);linspace(mean(dz_realtime),mean(dz_realtime),norbits)];

% Calculate Orbital Parameters

mu = 398600; %Gravitational Coefficient of Earth
maxorbitnumber=1;%Desired number of orbits until rendezvous

%Define Orbit Range and Step Size

tvalue=linspace(0.1*pi,maxorbitnumber*2*pi,10000);

%Find n value (Eq. 8.15 in Orbital Mechanics, Page 143, Conway)

%Define Orbit Parameters and Fill Proper Structure Fields
for i=1:norbits
    e=e(1);
    a = sqrt(r_(1)^2+r_(2)^2+r_(3)^2)/(1+e(1)); %Semi-Major Axis, km
    period_value=2*pi*sqrt(a^3/mu);
    deltaV_total(i).period=period_value;
    deltaV_total(i).dr_=dr_(:,i);
    deltaV_total(i).rstar_=rstar_(:,i);
    deltaV_total(i).n=sqrt(mu/(rstar_(1,i))^3);
    deltaV_total(i).Ra=sqrt(r_(1)^2+r_(2)^2+r_(3)^2);
    deltaV_total(i).Rp=sqrt(r_(1)^2+r_(2)^2+r_(3)^2);
end

%Find State Transition Matrix (Equation 8.26, page 149-150)
%Determine Initial DeltaV (Equations 8.34 and 8.36, page 151)
%Determine Final DeltaV (Equation 8.38, page 152)
%Find Total DeltaV (Combine Each)

%Note: Values above and below 0.25 and -.25 km/s are ignored


for j=1:norbits
    for i=1:length(tvalue)
	
	%Calculate Delta V
    
        tvaluetemp(j,i)=(tvalue(1,i)/2*pi)*deltaV_total(j).period;
        numberorbits(j,i)=(tvalue(1,i)/(2*pi));
        s=sin(deltaV_total(j).n*tvaluetemp(j,i)*(180/pi));
        c=cos(deltaV_total(j).n*tvaluetemp(j,i)*(180/pi));

        M=[(4-3*c) 0 0;6*(s-deltaV_total(j).n*tvalue(i)) 1 0;0 0 c];
        N=[s/deltaV_total(j).n (2/deltaV_total(j).n)*(1-c) 0; -(2/deltaV_total(j).n)*(1-c) (4*s-3*deltaV_total(j).n*tvaluetemp(j,i))/deltaV_total(j).n 0; 0 0 s/deltaV_total(j).n];
        S=[3*deltaV_total(j).n*s 0 0; -6*deltaV_total(j).n*(1-c) 0 0; 0 0 -deltaV_total(j).n*s];
        T=[c 2*s 0; -2*s 4*c-3 0; 0 0 c];

        deltavtempfinal=(T/N*M-S)*deltaV_total(:,j).dr_;
        deltavtempinitial=-inv(N)*M*deltaV_total(:,j).dr_-dv_initial(:,j);
        deltavtemp=deltavtempfinal+deltavtempinitial;
        if deltavtempfinal(1)>.25 || deltavtempinitial(1)>.25 || deltavtemp(1)>.25
            deltavx(j,i)=.25;
            deltavxi(j,i)=.25;
            deltavxf(j,i)=.25;
        elseif deltavtempfinal(1)<-.25 || deltavtempinitial(1)<-.25 || deltavtemp(1)<-.25
            deltavx(j,i)=-.25;    
            deltavxi(j,i)=-.25;
            deltavxf(j,i)=-.25;
        else
            deltavx(j,i)=deltavtemp(1);
            deltavxi(j,i)=deltavtempfinal(1);
            deltavxf(j,i)=deltavtempinitial(1);
        end
        
        if deltavtempfinal(2)>.25 || deltavtempinitial(2)>.25 || deltavtemp(2)>.25
            deltavy(j,i)=.25;
            deltavyi(j,i)=.25;
            deltavyf(j,i)=.25;
        elseif deltavtempfinal(2)<-.25 || deltavtempinitial(2)<-.25 || deltavtemp(2)<-.25
            deltavy(j,i)=-.25;    
            deltavyi(j,i)=-.25;
            deltavyf(j,i)=-.25;
        else
            deltavy(j,i)=deltavtemp(2);
            deltavyi(j,i)=deltavtempfinal(2);
            deltavyf(j,i)=deltavtempinitial(2);
        end
        
        if deltavtempfinal(3)>.25 || deltavtempinitial(3)>.25 || deltavtemp(3)>.25
            deltavz(j,i)=.25;
            deltavzi(j,i)=.25;
            deltavzf(j,i)=.25;
        elseif deltavtempfinal(3)<-.25 || deltavtempinitial(3)<-.25 || deltavtemp(3)<-.25
            deltavz(j,i)=-.25;    
            deltavzi(j,i)=-.25;
            deltavzf(j,i)=-.25;
        else
            deltavz(j,i)=deltavtemp(3);
            deltavzi(j,i)=deltavtempfinal(3);
            deltavzf(j,i)=deltavtempinitial(3);
        end
    end
    deltaV_total(j).numorbits=numberorbits(j,:);
    deltaV_total(j).t=tvaluetemp(j,:);
    deltaV_total(j).deltav_x=deltavx(j,:);
    deltaV_total(j).deltav_y=deltavy(j,:);
    deltaV_total(j).deltav_z=deltavz(j,:);
    deltaV_total(j).deltav_zi=deltavzi(j,:);
    deltaV_total(j).deltav_zf=deltavzf(j,:);
    deltaV_total(j).deltav_xi=deltavxi(j,:);
    deltaV_total(j).deltav_xf=deltavxf(j,:);
    deltaV_total(j).deltav_yi=deltavyi(j,:);
    deltaV_total(j).deltav_yf=deltavyf(j,:);
    deltaV_total(j).deltav_total=(deltavx(j,:).^2+deltavy(j,:).^2+deltavz(j,:).^2).^(.5);
end

%Create a Color Map

cc=jet(2*norbits);

%Plot Results

figure(1)

% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_x,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Delta V in Local X Direction (km/s)')
% title(['Total Delta Vx vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% figure(2)
% 
% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_xi,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Initial Delta V in Local X Direction (km/s)')
% title(['Initial Delta Vx vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% figure(3)
% 
% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_xf,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Final Delta V in Local X Direction (km/s)')
% title(['Final Delta Vx vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% figure(4)
% 
% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_y,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Delta V in Local Y Direction (km/s)')
% title(['Total Delta Vy vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% figure(5)
% 
% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_yi,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Initial Delta V in Local Y Direction (km/s)')
% title(['Initial Delta Vy vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% figure(6)
% 
% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_yf,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Final Delta V in Local Y Direction (km/s)')
% title(['Fianl Delta Vy vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% 
% figure(7)
% 
% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_z,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Delta V in Local Z Direction (km/s)')
% title(['Total Delta Vz vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% figure(8)
% 
% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_zi,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Initial Delta V in Local Z Direction (km/s)')
% title(['Initial Delta Vz vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% figure(9)
% 
% hold on
% for i=1:norbits
%     plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_zf,'color',cc(i,:))
%     legendinfo{i}=['Initial Separation (km): ' num2str(sqrt(3*dr_(1,i)^2))];
% end
% legend(legendinfo)
% xlabel('Number of Orbits Until Rendezvous')
% ylabel('Final Delta V in Local Z Direction (km/s)')
% title(['Final Delta Vz vs Number of Orbits Until Rendezvous, Starting Range between 0.001 km and ' num2str(sqrt(dr_(1,length(r_))^2+dr_(2,length(r_))^2+dr_(3,length(r_))^2)) ' km'])
% 
% figure(10)

hold on
for i=1:norbits
    plot(deltaV_total(i).numorbits,deltaV_total(i).deltav_total,'color',cc(i,:))
    legendinfo{i}=['Initial Relative Velocity (km/s): ' num2str(sqrt(3*dv_initial(1,i)^2))];
end
legend(legendinfo)
xlabel('Number of Orbits Until Rendezvous')
ylabel('Total Delta V (km/s)')
% title(['Total Delta V vs Number of Orbits Until Rendezvous, Starting Displacement of ' num2str(sqrt(dr_(1,norbits)^2+dr_(2,norbits)^2+dr_(3,norbits)^2)) ' km'])

    


















%=======================Rascal Senior Design===============================
%============================RCL-C-NAV2====================================
%=======================Author: Tom Moline=================================
%=====================Approval: Nate Richard===============================
%=====================Revision: - =========================================
%===============Last Edit Date: October 21, 2013===========================

%============================Summary=======================================
%This script serves as a preliminary method of solving lambert's problem,
%as described on page 78 of Orbital Mechanics, Prussing-Conway. The results
%of this function can be used as an input in the rendezvous_deltaV script
%to calculate the total deltaV requried to perform the desired orbital
%transfer.

%============================Inputs========================================
% r1_ = Inital Displacement Vector, [x;y;z] (km)
% r2_ = Final Displacemetn Vector, [x;y;z] (km)
% tf =Difference between time for target satellite and chasing satellite 
%     reach point 2 (km)
%
% Note: All variables with an underscore are vectors 
 

%===========================Outputs========================================
% Terminal velocity vectors at each point

%==========================Begin Code======================================

clear
clc

%Initialize Constants
mu=398600; %Earth's standard gravitational parameter,m^3/s^2
re=6371; %Earth's radius, km

%Initialize Postion and Velocity Vectors
r_ = [3212.59    4572    -3877.23]; % km
v_ = [-6.379    1.003    -4.106]; % km/s
dr_=[0 0 0]; %km
dv_=[0.001 0 0]; %km/s
orbitnumber=1;

[dr_realtime,dv_realtime,x,y,z,dx,dy,dz,dx_realtime,dy_realtime,dz_realtime,x_realtime,y_realtime,z_realtime,theta,T] = orbits_plot(r_,v_,dr_,dv_,orbitnumber);

t=1:round(orbitnumber*T(1));

figure(1)
plot(t./60,x_realtime,t./60,y_realtime,t./60,z_realtime)
xlabel('Time (min)')
ylabel('Relative Displacement (km)')
legend(['Relatvie Displacement in X Direction for: ', num2str(orbitnumber), ' orbits'],...
    ['Relatvie Displacement in Y Direction for: ', num2str(orbitnumber), ' orbits'],...
    ['Relatvie Displacement in Z Direction for: ', num2str(orbitnumber), ' orbits'])
    
figure(2)
plot(t./60,dx_realtime.*100000,t./60,dy_realtime.*100000,t./60,dz_realtime.*100000)
xlabel('Time (min)')
ylabel('Relative Velocity (cm/s)')
legend(['Relatvie Velocity in X Direction for: ', num2str(orbitnumber), ' orbits'],...
    ['Relatvie Velocity in Y Direction for: ', num2str(orbitnumber), ' orbits'],...
    ['Relatvie Velocity in Z Direction for: ', num2str(orbitnumber), ' orbits'])

figure(3)
plot(t./60,dr_realtime)
xlabel('Time (min)')
ylabel('Relative Displacement (km)')

figure(4)
plot(t./60,dv_realtime.*100000)
xlabel('Time (min)')
ylabel('Relative Velocity (cm/s)')


% N=0;
% 
% tf=95*60;
% 
% [minval,ind]=min(dr_realtime(:));
% [I,J]=ind2sub(size(dr_realtime),ind);
% 
% r1_=[x(1,J),y(1,J),z(1,J)];
% r2_=[x(2,J),y(2,J),z(2,J)];
% v1_=[dx(1,J),dy(1,J),dz(1,J)];
% v2_=[dx(2,J),dy(2,J),dz(2,J)];
% 
% r1=sqrt(dot(r1_,r1_));
% r2=sqrt(dot(r2_,r2_));
% 
% s=(sqrt(x(1,J)^2+y(1,J)^2+z(1,J)^2)+sqrt(x(2,J)^2+y(2,J)^2+z(2,J)^2)+dr_realtime(J))/2;
% 
% c=dr_realtime(J);
% m=r1+r2+c;
% n=r1+r2-c;
% 
% alpha=dot([0 0 1],cross(r1_,r2_));
% psi_0=acos(dot(r1_,r2_)/(r1*r2));
% 
% psi=2*pi-psi_0;
% 
% sigma=sqrt(((4*r1*r2)/m^2)*(cos(psi/2))^2);
% 
% if sigma^2 ~= 0
%     if sigma >0
%         sign_sigma=1;
%     else
%         sign_sigma=2;
%     end
% else
%     sign_sigma=1;
% end
% 
% tau=4*tf*sqrt(mu/m^3);
% taup=(2/3)*(1-sigma^3);
% 
% Nmax=round(tau/pi);
% 
% tau_me=N*pi+acos(sigma)+sigma*sqrt(1-sigma^2);
% 
% if tau < tau_me
%     x1=.5;
% elseif tau==tau_me
%     x1=0;
% else
%     x1=-.5;
% end
% 
% if sign_sigma==1
%     y1=sqrt(1-sigma^2*(1-x1^2));
% else
%     y1=-sqrt(1-sigma^2*(1-x1^2));
% end
% 
% phi_x0=acot(x1/sqrt(1-x1^2))-(1/(3*x1))*(2+x1^2)*sqrt(1-x1^2);
% phi_y0=acot(y1/sqrt(1-y1^2))-(1/(3*y1))*(2+y1^2)*sqrt(1-y1^2);
% 
% F_0=phi_x0+phi_y0+N*pi-tau;
% 
% x2=x1-((n*F_0)/((((sigma^2*x1)/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 1)^(1/2))/(-sigma^2*(x1^2 - 1))^(3/2))/(((x1^2 - 1)*sigma^2 + 1)/(sigma^2*(x1^2 - 1)) - 1) + (x1^2/(1 - x1^2)^(3/2) + 1/(1 - x1^2)^(1/2))/(x1^2/(x1^2 - 1) - 1) + (x1^2 + 2)/(2*(1 - x1^2)^(1/2)) - (1 - x1^2)^(1/2) + ((1 - x1^2)^(1/2)*(x1^2 + 2))/(2*x1^2) - (sigma^2*x1*(-sigma^2*(x1^2 - 1))^(1/2))/((x1^2 - 1)*sigma^2 + 1)^(1/2) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 3))/(2*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 3))/(2*((x1^2 - 1)*sigma^2 + 1)^(3/2)))...
%     +(((sigma^2*x1)/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 1)^(1/2))/(-sigma^2*(x1^2 - 1))^(3/2))/(((x1^2 - 1)*sigma^2 + 1)/(sigma^2*(x1^2 - 1)) - 1) + (x1^2/(1 - x1^2)^(3/2) + 1/(1 - x1^2)^(1/2))/(x1^2/(x1^2 - 1) - 1) + (x1^2 + 2)/(2*(1 - x1^2)^(1/2)) - (1 - x1^2)^(1/2) + ((1 - x1^2)^(1/2)*(x1^2 + 2))/(2*x1^2) - (sigma^2*x1*(-sigma^2*(x1^2 - 1))^(1/2))/((x1^2 - 1)*sigma^2 + 1)^(1/2) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 3))/(2*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 3))/(2*((x1^2 - 1)*sigma^2 + 1)^(3/2)))...
%     /abs(((sigma^2*x1)/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 1)^(1/2))/(-sigma^2*(x1^2 - 1))^(3/2))/(((x1^2 - 1)*sigma^2 + 1)/(sigma^2*(x1^2 - 1)) - 1) + (x1^2/(1 - x1^2)^(3/2) + 1/(1 - x1^2)^(1/2))/(x1^2/(x1^2 - 1) - 1) + (x1^2 + 2)/(2*(1 - x1^2)^(1/2)) - (1 - x1^2)^(1/2) + ((1 - x1^2)^(1/2)*(x1^2 + 2))/(2*x1^2) - (sigma^2*x1*(-sigma^2*(x1^2 - 1))^(1/2))/((x1^2 - 1)*sigma^2 + 1)^(1/2) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 3))/(2*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 3))/(2*((x1^2 - 1)*sigma^2 + 1)^(3/2)))...
%     *sqrt((n-1^2)*(((sigma^2*x1)/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 1)^(1/2))/(-sigma^2*(x1^2 - 1))^(3/2))/(((x1^2 - 1)*sigma^2 + 1)/(sigma^2*(x1^2 - 1)) - 1) + (x1^2/(1 - x1^2)^(3/2) + 1/(1 - x1^2)^(1/2))/(x1^2/(x1^2 - 1) - 1) + (x1^2 + 2)/(2*(1 - x1^2)^(1/2)) - (1 - x1^2)^(1/2) + ((1 - x1^2)^(1/2)*(x1^2 + 2))/(2*x1^2) - (sigma^2*x1*(-sigma^2*(x1^2 - 1))^(1/2))/((x1^2 - 1)*sigma^2 + 1)^(1/2) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 3))/(2*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 3))/(2*((x1^2 - 1)*sigma^2 + 1)^(3/2)))^2-...
%     n*(n-1)*F_0*(1 - x1^2)^(1/2)/x1 + (sigma^2/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*((x1^2 - 1)*sigma^2 + 1)^(1/2))/(-sigma^2*(x1^2 - 1))^(3/2) - (sigma^4*x1^2)/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(3/2)) + (2*sigma^4*x1^2)/((-sigma^2*(x1^2 - 1))^(3/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (3*sigma^4*x1^2*((x1^2 - 1)*sigma^2 + 1)^(1/2))/(-sigma^2*(x1^2 - 1))^(5/2))/(((x1^2 - 1)*sigma^2 + 1)/(sigma^2*(x1^2 - 1)) - 1) + ((3*x1^3)/(1 - x1^2)^(5/2) + (3*x1)/(1 - x1^2)^(3/2))/(x1^2/(x1^2 - 1) - 1) + (2*x1)/(1 - x1^2)^(1/2) + (x1*(x1^2 + 2))/(2*(1 - x1^2)^(3/2)) - ((x1^2/(1 - x1^2)^(3/2) + 1/(1 - x1^2)^(1/2))*((2*x1)/(x1^2 - 1) - (2*x1^3)/(x1^2 - 1)^2))/(x1^2/(x1^2 - 1) - 1)^2 - (((sigma^2*x1)/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*x1*((x1^2 - 1)*sigma^2 + 1)^(1/2))/(-sigma^2*(x1^2 - 1))^(3/2))*((2*x1)/(x1^2 - 1) - (2*x1*((x1^2 - 1)*sigma^2 + 1))/(sigma^2*(x1^2 - 1)^2)))/(((x1^2 - 1)*sigma^2 + 1)/(sigma^2*(x1^2 - 1)) - 1)^2 - (sigma^2*(-sigma^2*(x1^2 - 1))^(1/2))/((x1^2 - 1)*sigma^2 + 1)^(1/2) - (x1^2 + 2)/(2*x1*(1 - x1^2)^(1/2)) - ((1 - x1^2)^(1/2)*(x1^2 + 2))/x1^3 + (sigma^2*((x1^2 - 1)*sigma^2 + 3))/(2*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (sigma^2*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 3))/(2*((x1^2 - 1)*sigma^2 + 1)^(3/2)) + (2*sigma^4*x1^2)/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) + (2*sigma^4*x1^2*(-sigma^2*(x1^2 - 1))^(1/2))/((x1^2 - 1)*sigma^2 + 1)^(3/2) - (sigma^4*x1^2*((x1^2 - 1)*sigma^2 + 3))/((-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 1)^(3/2)) + (sigma^4*x1^2*((x1^2 - 1)*sigma^2 + 3))/(2*(-sigma^2*(x1^2 - 1))^(3/2)*((x1^2 - 1)*sigma^2 + 1)^(1/2)) - (3*sigma^4*x1^2*(-sigma^2*(x1^2 - 1))^(1/2)*((x1^2 - 1)*sigma^2 + 3))/(2*((x1^2 - 1)*sigma^2 + 1)^(5/2)))))l;






        


















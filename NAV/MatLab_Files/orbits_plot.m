%=======================Rascal Senior Design===============================
%============================RCL-C-NAV3====================================
%=======================Author: Tom Moline=================================
%=====================Approval: Nate Richard===============================
%=====================Revision: - =========================================
%===============Last Edit Date: October 21, 2013===========================

%============================Summary=======================================
%This script plots the orbits of two spacecraft over a specified period of
%time, provided an inital relative velocity vector and relative position.

%============================Inputs========================================
% rd_ = Inital Relative Displacement Vector, [x;y;z] (km)
% vd_ = Initial Relative Velocity Vector, [vx;vy;vz] (km/s)
% r_=Inital Displacement Vector Relative to Earth Focal Point [x;y;z](km) 
% v_=Inital Velocity Vector Relative to Earth Focal Point [vx;vy;vz] (km/s)
% t = Time over Which to Plot Results
%
% Note: All variables with an underscore are vectors 
 

%===========================Outputs========================================
% Plot of orbits over a specified time

%==========================Begin Code======================================

clear
clc

%========================Define Constants==================================

mew = 398600; %Gravity Constant for Earth
re=6371; % Radius of Earth, km

%=========================Define Inital Conditions=========================

r_ = [6863.700    5877.350    -3347.670]; % km
v_ = [-5.1039    3.2546    4.68185]; % km/s
dr_=[0 0 0]; %km
dv_=[0.0005 0 0]; %km/s

for i=1:3
    r_(2,i)=r_(1,i)+dr_(1,i);
    v_(2,i)=v_(1,i)+dv_(1,i);
end

%======================Find Period of Orbit================================
for i=1:2
    r(i) = sqrt( r_(i,1)^2+r_(i,2)^2+r_(i,3)^2);
    a(i)= -mew * r(i)/((v_(i,1)^2+v_(i,2)^2+v_(i,3)^2)*r(i)-2*mew);
    T(i)=(2*pi)/sqrt(mew)*a(i)^(3/2);
    t(i)=T(i)/1000;
end

%===========Find Magnitude of Positon and Velocity Vectors=================
% r = sqrt( r_(1)^2+r_(2)^2+r_(3)^2);
% v=sqrt(dot(v_,v_));
% Vr = dot(v_, r_)/r;
% 
% %===================Find Specific Momentum Vector==========================
% h_ = cross(r_,v_);
% h = sqrt(dot(h_,h_));
% 
% %========================Find Eccentricity=================================
% 
% e_ = (1/mew)*((v^2-(mew/r))*r_-r*Vr*v_);
% e=sqrt(dot(e_,e_));
% 
% 
% %===========================Find Semi-Major Axis===========================
% a= -mew * r/(dot(v_,v_) * r - 2*mew);
% N_ = cross([0 0 1], h_);
% N = sqrt(dot(N_,N_));
% 
% %========================Find Argument of Perigee==========================
% w = acos(dot(N_,e_)/(N*e));
% if(e_(3) < 0)
%     w=2*pi() - w;
% end
% N = sqrt(dot(N_,N_));
% 
% %=========================Find True Anomaly================================
% theta= acos(dot(e_,r_)/(e*r));
% if (theta < 0)
%     theta = 2*pi-theta;
% end
% 
% %===========================Find Inclination===============================
% i = acos( h_(3) / h);
% 
% %===============Find Longitude of Ascending Node===========================
% if( N_(2) < 0)
%     omega = 2*pi() - acos(N_(1)/N);
% else
%     omega = acos(N_(1)/N);
% end
% 
% %========================Find Mean Anomaly=================================
% Mnot = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2))-(e*sqrt(1-e^2)*sin(theta))/(1+e*cos(theta));
% Enot = acos((e+cos(theta))/(1+e*cos(theta)));
% E=Enot;
% M=Mnot;
% for k=1:5*T
% %========================Find Eccentric Anomaly============================
%     E_ratio=1;
%     M = M + sqrt(mew/a^3)*(t);
%     if M<pi
%         E=M+e/2;
%     else
%         E=M-e/2;
%     end
%     
%     while E_ratio>10^-8
%         E_new=E-e*sin(E)-M;
%         E_prime=1-e*cos(E);
%         E_ratio=E_new/E_prime;
%         E=E-E_ratio;
%     end
%     theta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
%     rp_=(h^2/mew)*(1/(1+e*cos(theta)))*[cos(theta);sin(theta);0];
%     vp_=(mew/h)*[-sin(theta);e+cos(theta);0];
%     
% %=========================Create Perifocal Transformation Matrix===========
%     QXx_=[-sin(omega)*cos(i)*sin(w)+cos(omega)*cos(w), cos(omega)*cos(i)*sin(w)+sin(omega)*cos(w), sin(i)*sin(w);
%         -sin(omega)*cos(i)*cos(w)-cos(omega)*sin(w), cos(omega)*cos(i)*cos(w)-sin(omega)*sin(w), sin(i)*cos(w);
%         sin(omega)*sin(i), -cos(omega)*sin(i), cos(i)];
%     
% %===========================Create Geocentric Transformation Matrix========
%     QxX_=[-sin(omega)*cos(i)*sin(w)+cos(omega)*cos(w), -sin(omega)*cos(i)*cos(w)-cos(omega)*sin(w), sin(omega)*sin(i);
%         cos(omega)*cos(i)*sin(w)+sin(omega)*cos(w), cos(omega)*cos(i)*cos(w)-sin(omega)*sin(w), -cos(omega)*sin(i);
%         sin(w)*sin(i), cos(w)*sin(i), cos(i)];
%     
% %==============================Define Perifocal Quantities=================
%     r_temp=QxX_*rp_;
%     v_temp=QxX_*vp_;
%     
%     r_=r_temp.';
%     v_=v_temp.';
%     
%     x(k)=r_(1);
%     y(k)=r_(2);
%     z(k)=r_(3);
%     
% end
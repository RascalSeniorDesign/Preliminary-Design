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

r_ = [3212.59    4572    -3877.23]; % km
v_ = [-6.379    1.003    -4.106]; % km/s
dr_=[0 0 0]; %km
dv_=[0.0005 .0005 .0005]; %km/s

for i=1:3
    r_(2,i)=r_(1,i)+dr_(1,i);
    v_(2,i)=v_(1,i)+dv_(1,i);
end

%======================Find Period of Orbit================================
for i=1:2
    r(i) = sqrt( r_(i,1)^2+r_(i,2)^2+r_(i,3)^2);
    a(i)= r(i)/(2-((r(i)*(v_(i,1)^2+v_(i,2)^2+v_(i,3)^2))/mew));
    T(i)=(2*pi)/sqrt(mew)*a(i)^(3/2);
    t(i)=T(i)/1000;
end

for i=1:2
%===================Find Specific Momentum Vector==========================
    tempr_=[r_(i,1) r_(i,2) r_(i,3)];
    tempv_=[v_(i,1) v_(i,2) v_(i,3)];
    v(i)=sqrt(dot(tempv_,tempv_));
    htemp_ = cross(tempr_,tempv_);
    h(i) = sqrt(dot(htemp_,htemp_));
%========================Find Eccentricity=================================
    tempe_ = (1/mew)*((v(i)^2-(mew/r(i)))*tempr_-dot(tempr_,tempv_)*tempv_);
    e(i)=sqrt(dot(tempe_,tempe_));
%===========================Find Semi-Major Axis===========================
    tempN_ = cross([0 0 1], htemp_);
    N(i) = sqrt(dot(tempN_,tempN_));

%========================Find Argument of Perigee==========================
    if(tempe_(3) < 0)
        w(i)=2*pi() - w;
    else
        w(i) = acos(dot(tempN_,tempe_)/(N(i)*e(i)));
    end

%=========================Find True Anomaly================================
    theta(i)= acos(dot(tempe_,tempr_)/(e(i)*r(i)));
    if (theta < 0)
        theta(i) = 2*pi-theta;
    else
        theta(i)=theta(i);
    end
%===========================Find Inclination===============================
    inc(i) = acos( htemp_(3) / h(i));
%===============Find Longitude of Ascending Node===========================
    if( tempN_(2) < 0)
        omega(i) = 2*pi - acos(tempN_(1)/N(i));
    else
        omega(i) = acos(tempN_(1)/N(i));
    end

%========================Find Mean Anomaly=================================
    Mnot(i) = 2*atan(sqrt((1-e(i))/(1+e(i)))*tan(theta(i)/2))-(e(i)*sqrt(1-e(i)^2)*sin(theta(i)))/(1+e(i)*cos(theta(i)));
    Enot(i) = acos((e(i)+cos(theta(i)))/(1+e(i)*cos(theta(i))));
    E=Enot(i);
    M=Mnot(i);
    
%=========================Create Perifocal Transformation Matrix===========
    QXx_=[-sin(omega(i))*cos(inc(i))*sin(w(i))+cos(omega(i))*cos(w(i)), cos(omega(i))*cos(inc(i))*sin(w(i))+sin(omega(i))*cos(w(i)), sin(inc(i))*sin(w(i));...
            -sin(omega(i))*cos(inc(i))*cos(w(i))-cos(omega(i))*sin(w(i)), cos(omega(i))*cos(inc(i))*cos(w(i))-sin(omega(i))*sin(w(i)), sin(inc(i))*cos(w(i));...
            sin(omega(i))*sin(inc(i)), -cos(omega(i))*sin(inc(i)), cos(inc(i))];
    
%===========================Create Geocentric Transformation Matrix========
    QxX_=[-sin(omega(i))*cos(inc(i))*sin(w(i))+cos(omega(i))*cos(w(i)), -sin(omega(i))*cos(inc(i))*cos(w(i))-cos(omega(i))*sin(w(i)), sin(omega(i))*sin(inc(i));...
            cos(omega(i))*cos(inc(i))*sin(w(i))+sin(omega(i))*cos(w(i)), cos(omega(i))*cos(inc(i))*cos(w(i))-sin(omega(i))*sin(w(i)), -cos(omega(i))*sin(inc(i));...
            sin(w(i))*sin(inc(i)), cos(w(i))*sin(inc(i)), cos(inc(i))];
    

    for k=1:T(1)
%========================Find Eccentric Anomaly============================
        E_ratio=1;
        M = M + sqrt(mew/a(i)^3)*(t(i));
        if M<pi
            E=M+e/2;
        else
            E=M-e/2; 
        end
    
        while E_ratio>10^-8
            E_new=E-e(i)*sin(E)-M;
            E_prime=1-e(i)*cos(E);
            E_ratio=E_new/E_prime;
            E=E-E_ratio;
        end
        thetatemp=2*atan(sqrt((1+e(i))/(1-e(i)))*tan(E(1)/2));
        rp_=(h(i)^2/mew)*(1/(1+e(i)*cos(thetatemp)))*[cos(thetatemp);sin(thetatemp);0];
        vp_=(mew/h(i))*[-sin(thetatemp);e(i)+cos(thetatemp);0];

%==============================Define Perifocal Quantities=================
        r_temp=QxX_*rp_;
        v_temp=QxX_*vp_;
    
        rtemp_=r_temp.';
        vtemp_=v_temp.';
    
        x(i,k)=rtemp_(1);
        y(i,k)=rtemp_(2);
        z(i,k)=rtemp_(3);
        
        dx(i,k)=vtemp_(1);
        dy(i,k)=vtemp_(2);
        dz(i,k)=vtemp_(3);
    end
end

for k=1:T(1)
    dr_realtime(k)=((x(1,k)-x(2,k))^2+(y(1,k)-y(2,k))^2+(z(1,k)-z(2,k))^2)^(1/2);
    dv_realtime(k)=((dx(1,k)-dx(2,k))^2+(dy(1,k)-dy(2,k))^2+(dz(1,k)-dz(2,k))^2)^(1/2);
    dx_realtime(k)=dx(1,k)-dx(2,k);
    dy_realtime(k)=dy(1,k)-dy(2,k);
    dz_realtime(k)=dz(1,k)-dz(2,k);
end

%cc=jet(2);
satellite=['Target';'Chaser'];

tstep=1:T(1);


figure(1)
plot(tstep./60,dr_realtime)

figure(2)
plot(tstep./60,dv_realtime.*100000)

figure(3)
plot(tstep./60,dx_realtime.*100000,tstep./60,dy_realtime.*100000,tstep./60,dz_realtime.*100000)


% % for j=1:tstep:T(1)/2
% %         b=round(1+j*t/tstep);
% %         set(Target,'matrix',makehgtform('translate',x(1,b),y(1,b),z(1,b)));
% %         set(Chaser,'matrix',makehgtform('translate',x(2,b),y(2,b),z(2,b)));
% %         line(x(1,1:b),y(1,1:b),z(1,1:b),'color','b')
% %         line(x(2,1:b),y(2,1:b),z(2,1:b),'color','r')
% %         drawnow
% %     legendinfo{i}=[satellite(i,:) ' Satellite Orbit for Initial Relative Velocity of ' num2str(dv_(1)*100000) 'cm/s in x Direction'];
% % end

% legend(legendinfo)
% xlabel('X Position (km)')
% ylabel('Y Position (km)')
% zlabel('Z Position (km)')



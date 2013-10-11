%=======================Rascal Senior Design===============================
%============================RCL-C-NAV1====================================
%=======================Author: Tom Moline=================================
%=====================Approval: Nate Richard===============================
%=====================Revision: - =========================================
%===============Last Edit Date: October 9, 2013============================

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

Rp = input('Please input the periapsis of a selected orbit: ');
Ra = input('Please input the apoapsis of a selected orbit: ');

% Calculate Orbital Parameters
e = (Ra-Rp)/(Ra+Rp); %Eccentricity of Orbit
a = Ra/(1+e); %Semi-Major Axis, km
mu = 398600; %Standard Gravitational Parameter of Earth, m^3/s^2
period=2*pi*sqrt(a^3/mu);

% Ask for Initial Locations of Satellites
rstar_ = input('Please input the inital displacement vector of chase satellite: ');
r_ = input('Please input the initial displacement vector of target satellite: ');

% Calculate Local Gravity Gradient Matrix

grav_=(mu/rstar_(1)^3)*[2 0 0 ; 0 -1 0; 0 0 -1];
n=sqrt(mu/rstar_(1)^3);

%Indicate desired orbital time step for final rendezvous

t=linspace(.5,12*pi,50);

for i=1:length(t)
	
	% Calculate Delta V

	s=sin(n*t(i)*(180/pi));
	c=cos(n*t(i)*(180/pi));

	M=[4-3*c 0 0; 6*(s-n*t) 1 0; 0 0 c];
	N=[s/n (2/n)*(1-c) 0; -(2/n)*(1-c) (4*s-3*n*t)/n 0; 0 0 s/n];
	S=[3*n*s 0 0; -6*n*(1-c) 0 0; 0 0 -n*s];
	T=[c 2*s 0; -2*s 4*c-3 0; 0 0 c];

	deltav_x(i)=(T*inv(N)*M-S)*r_(1);
	deltav_y(i)=(T*inv(N)*M-S)*r_(2);

end

plot(t,deltav_x,t,deltav_y)
xlabel('Time (Period Step Size)')
ylabel('Delta V (km/s)')
legend('Delta V: X Direction','Delta V: Y Direction')















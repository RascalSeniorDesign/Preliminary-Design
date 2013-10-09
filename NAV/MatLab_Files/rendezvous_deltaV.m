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
% delta_r0 =Inital Relative Displacement Vector b/w Sundance and Butch(km)
%    NOTE: This is assumed to be a 2-D Vector input [x, y]

%===========================Outputs========================================
% Total Delta V to go from intial positon to one within 10 m of Sundance

%==========================Begin Code======================================

Rp = input('Please input the periapsis of a selected orbit: ');
Ra = input('Please input the apoapsis of a selected orbit: ');

% Calculate Orbital Parameters
e = (Ra-Rp)/(Ra+Rp); %Eccentricity of Orbit
a = Ra/(1+e); %Semi-Major Axis, km
mu = 398600; %Standard Gravitational Parameter of Earth, m^3/s^2



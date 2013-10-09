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
%    NOTE: This is assumed to be a 2-D Vector input [x, y]
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

% Ask for Initial Locations of Satellites
rstar_ = input('Please input the inital displacement vector of chase satellite: ');
r_ = input('Please input the initial displacement vector of target satellite: ');











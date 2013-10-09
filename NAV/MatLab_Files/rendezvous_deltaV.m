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

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
% r1_=Inital Displacement Vector Relative to Earth Focal Point [x;y;z](km) 
% t = Time over Which to Plot Results
%
% Note: All variables with an underscore are vectors 
 

%===========================Outputs========================================
% Plot of orbits over a specified time

%==========================Begin Code======================================
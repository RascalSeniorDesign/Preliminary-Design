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

%Initialize Constants
mu=398600; %Earth's standard gravitational parameter,m^3/s^2
re=6371; %Earth's radius, km

%Initialize Postion Vectors and Final Time
r1_=[re+500;0;0];
r2_=[re+500.1;.1;0];


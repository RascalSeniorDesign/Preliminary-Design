%=======================Rascal Senior Design===============================
%============================RCL-C-NAV2====================================
%=======================Author: Tom Moline=================================
%=====================Approval: Nate Richard===============================
%=====================Revision: - =========================================
%===============Last Edit Date: October 13, 2013===========================

%============================Summary=======================================
%This script serves as a preliminary method of solving lambert's problem,
%as described on page 78 of Orbital Mechanics, Prussing-Conway. The results
%of this function can be used as an input in the rendezvous_deltaV script
%to calculate the total deltaV requried to perform the desired orbital
%transfer.

%============================Inputs========================================
% R1_ = Inital Displacement Vector, [x;y;z] (km)
% R2_ = Final Displacemetn Vector, [x;y;z] (km)
% tf =Difference between time for target satellite and chasing satellite 
%     reach point 2 (km)
%
% Note: All variables with an underscore are vectors 
 

%===========================Outputs========================================
% Terminal velocity vectors at each point

%==========================Begin Code======================================

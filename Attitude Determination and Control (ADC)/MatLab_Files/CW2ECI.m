function [Xint_] = CW2ECI(Xtgt_,x,y,z,xdot,ydot,zdot)
%The CW2ECI function, when provided an inertial target position and
%relative interceptors positions/velocities and finds the inertial position
%and velocity of the interceptor.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      Xtgt_     Inertial Target Pos/Velocity 2x3 Matrix      km and km/s
%      x         Initial Relative X Position   Scalar             km
%      y         Initial Relative Y Position   Scalar             km
%      z         Initial Relative Z Position   Scalar             km
%      xdot      Initial Relative X Velocity   Scalar            km/s
%      ydot      Initial Relative Y Velocity   Scalar            km/s
%      zdot      Initial Relative Z Velocity   Scalar            km/s
%      Xint_     Inertial Int Pos/Velocity   2x3 Matrix       km and km/s
%==========================================================================
%Initial Release, CW2ECI.m, Tom Moline, 2/5/2014


%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
rtgt_=Xtgt_(:,1);
vtgt_=Xtgt_(:,2);
rtgt=sqrt(sum(abs(rtgt_).^2));
R=rtgt_./rtgt; %R unit vector in local CW frame, target centered
W=cross(rtgt_,vtgt_)./sqrt(sum(abs(cross(rtgt_,vtgt_)).^2)); %W unit vector
S=cross(W,R); %S unit vector
vhill_=[xdot;ydot;zdot];
rhill_=[x;y;z];
RSW=[R S W];


XtgtRSW_=RSW'*Xtgt_;
rtgtRSW_=XtgtRSW_(:,1);
vtgtRSW_=XtgtRSW_(:,2);
rtgtRSW=sqrt(sum(abs(rtgtRSW_).^2));
rinttemp1_=[rtgtRSW_(1)+x;rtgtRSW_(2);rtgtRSW_(3)];
rvtgtCross_=cross(rtgt_,vtgt_);
vintRSW_=vtgtRSW_+vhill_+cross((rvtgtCross_/(rtgt^2)),rhill_);
deltaAlphaY=y/rtgtRSW;

%Define rotation functions to be used for finding realtive interceptor pos.
rot2= @(theta) [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
rot3= @(beta)  [cos(beta) sin(beta) 0; -sin(beta) cos(beta) 0; 0 0 1];

rinttemp2_=rot3(-deltaAlphaY)*rinttemp1_;
deltaAlphaZ=z/rtgtRSW;

rintRSW_=rot2(deltaAlphaZ)*rinttemp2_;
 
 XintRSW_=[rintRSW_ vintRSW_];
 
 Xint_=RSW*XintRSW_;
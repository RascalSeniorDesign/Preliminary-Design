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
R=[1;0;0];
S=[0;1;0];
W=[0;0;1];
rtgt_=Xtgt_(1,:);
vtgt_=Xtgt_(2,:);
vhill_=xdot+ydot+zdot;
rhill_=x+y+z;
rtgt=sqrt(sum(abs(rgt_).^2));

XtgtRSW_=[R S W]'*Xtgt_;
rtgtRSW_=XtgtRSW_(1,:);
vtgtRSW_=XtgtRSW_(2,:);
rinttemp1_=rtgtRSW_+x;
vintRSW_=vtgtRSW_+vhill_+cross((cross(rtgt_,vtgt_)/(sum(abs(rtgt_).^2))),rhill_);
deltaAlphaY=y/rtgt;

rinttemp2_=[cos(-deltaAlphaY) sin(-deltaAlphaY) 0; -sin(-deltaAlphaY)...
    cos(-deltaAlphaY) 0; 0 0 1]*rinttemp1_;
deltaAlphaZ=z/rtgt;

rintRSW_=[cos(dletaAlphaZ) 0 -sin(deltaAlphaZ); 0 1 0; sind(deltaAlphaX)...
     0 cos(deltaAlphaZ)]*rinttemp2_;
 
 XintRSW_=[rintRSW_;vintRSW_];
 
 Xint_=[R S W]*XintRSW_;
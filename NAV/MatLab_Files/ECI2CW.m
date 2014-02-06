function [x,y,z,xdot,ydot,zdot] = ECI2CW(Xtgt_,Xint_)
%The ECI2CW function, when provided an inertial target and interceptor
%positions and velocities outputs the relative distance and velocity of a
%interceptor realtive to a target satellite in the target satellites local
%RSW frame.
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
%Initial Release, ECI2CW.m, Tom Moline, 2/5/2014


%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
rtgt_=Xtgt_(:,1);
vtgt_=Xtgt_(:,2);
rtgt=sqrt(sum(abs(rtgt_).^2));

R=rtgt_'/rtgt;
W=cross(rtgt_',vtgt_')/sqrt(sum(abs(cross(rtgt_',vtgt_')).^2));
S=cross(W,R);
RSW=[R;S;W];

XtgtRSW_=RSW'*Xtgt_;
XintRSW_=RSW'*Xint_;

rtgtRSW_=XtgtRSW_(:,1);
vtgtRSW_=XtgtRSW_(:,2);
rintRSW_=XintRSW_(:,1);
vintRSW_=XintRSW_(:,2);
rtgtRSW=sqrt(sum(abs(rtgtRSW_').^2));

deltaAlphaz=atan(rintRSW_(3,1)/rtgt);
rinttem_=[cos(-deltaAlphaz) 0 -sin(-deltaAlphaz); 0 1 0; sind(-deltaAlphaz)...
     0 cos(-deltaAlphaz)]*rintRSW_;
 deltaAlphaY=atan(rintRSW_(2,1)/rtgt);
 
 rinttemp_=[cos(deltaAlphaY) sin(deltaAlphaY) 0; -sin(deltaAlphaY)...
    cos(deltaAlphaY) 0; 0 0 1]*rinttem_;

rinthill_=[rinttemp_(1,1)-rtgtRSW(1,1); deltaAlphaY*rtgt; deltaAlphaz*rtgt];
rinthill_=rinthill_';

vinthill_=vintRSW_'-vtgtRSW_'-cross((cross(rtgtRSW_',vtgtRSW_')/rtgtRSW^2),rinthill_');

x=rinthill_(1,1);
y=rinthill_(1,2);
z=rinthill_(1,3);

xdot=vinthill_(1,1);
ydot=vinthill_(1,2);
zdot=vinthill_(1,3);


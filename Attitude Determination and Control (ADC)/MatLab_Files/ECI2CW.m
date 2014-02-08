function [x,y,z,xdot,ydot,zdot] = ECI2CW(Xtgt_,Xint_)
%The ECI to CW function converts inetial earth centric positons and
%velocities of a target and interceptor spacecraft to RSW coordinates
%relative to the target spacecraft
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      Xtgt_     Inertial Target Pos/Velocity 3x2 Matrix      km and km/s
%      x         Initial Relative X Position   Scalar             km
%      y         Initial Relative Y Position   Scalar             km
%      z         Initial Relative Z Position   Scalar             km
%      xdot      Initial Relative X Velocity   Scalar            km/s
%      ydot      Initial Relative Y Velocity   Scalar            km/s
%      zdot      Initial Relative Z Velocity   Scalar            km/s
%      Xint_     Inertial Int Pos/Velocity   3x2 Matrix       km and km/s
%==========================================================================
%Initial Release, ECI2CW.m, Tom Moline, 2/7/2014


%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
rtgt_=Xtgt_(:,1); %Target ECI position, km
vtgt_=Xtgt_(:,2); %Target ECI Velcity, km
rtgt=sqrt(sum(abs(rtgt_).^2)); %Magnitude of Target position, km

R=rtgt_./rtgt; %R unit vector in local CW frame, target centered
W=cross(rtgt_,vtgt_)./sqrt(sum(abs(cross(rtgt_,vtgt_)).^2)); %W unit vector
S=cross(W,R); %S unit vector

RSW=[R;S;W]; %RSW transformation matrix

%==========================================================================
%                    Begin Coordinate Conversion
%==========================================================================
XtgtRSW_=RSW'*Xtgt_; %Transform from ECI to RSW coordinates
XintRSW_=RSW'*Xint_;

%==========================================================================
%                    Begin Finding Relative Position
%==========================================================================
deltaZ=atan(XintRSW_(3,1)/rtgt); %Set up to find 

%Define rotation functions to be used for finding realtive interceptor pos.
rot2= @(theta) [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
rot3= @(beta)  [cos(beta) sin(beta) 0; -sin(beta) cos(beta) 0; 0 0 1];

%Define intial rotation
rintTem_=rot2(-deltaZ)*XintRSW_(:,1);

%Define second rotation
deltaY=atan(XintRSW_(2,1)/rtgt);
rintTemp_=rot3(deltaY)*rintTem_;

%==========================================================================
%                  Find Relative Position and Velocity
%==========================================================================
rintHill_=[rintTemp_(1,1)-XtgtRSW_(1,1); deltaY*rtgt; deltaZ*rtgt];
vintHill_=XintRSW_(:,2)-XtgtRSW_(:,2)-cross(cross(XtgtRSW_(:,1),XtgtRSW_(:,2))...
    /sum(abs(XtgtRSW_(:,1)).^2),rintHill_);

x=rintHill_(1,1); %Relative transverse track position
y=rintHill_(2,1); %Relative in track position
z=rintHill_(3,1); %Relative out of plane position
xdot=vintHill_(1,1); %Relative velocity components
ydot=vintHill_(2,1);
zdot=vintHill_(3,1);
end

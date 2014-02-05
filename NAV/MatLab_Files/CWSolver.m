function [x,y,z,xdot,ydot,zdot] = CW2ECI(x0,y0,z0,x0dot,y0dot,z0dot,rtgt,tf)
%The CWSolver function takes in the intial realtive position and velocity
%of an interceptor satellite relative to a target spacecraft that is in a
%circular orbit defined by omegatgt. It outputs the relative positoin and 
%velocity of the interceprotr with respect to the target after some time
%tf. This function assumes near-circular orbits and short time intervals
%for calculations.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      x0        Initial Relative X Positon    Scalar             km
%      y0        Initial Relative Y Position   Scalar             km
%      z0        Initial Relative Z Position   Scalar             km
%      x0dot     Initial Relative X Velocity   Scalar            km/s
%      y0dot     Initial Relative Y Velocity   Scalar            km/s
%      z0dot     Initial Relative Z Velocity   Scalar            km/s
%      tf        Time Final Positoin Time      Scalar             min
%==========================================================================
%Initial Release, CWSolver.m, Tom Moline, 1/29/2014


%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
% x0=x0/6378.137;
% y0=y0/6378.137;
% z0=z0/6378.137;
% rtgt=rtgt/6378.137;
% x0dot=x0dot/7.9053838;
% y0dot=y0dot/7.9053838;
% z0dot=z0dot/7.9053838;
% tf=tf/13.446849;

tf=tf*60;

omegatgt=sqrt(398600.5/rtgt^3);

%==========================================================================
%                       Update to tf Position
%==========================================================================
x=(x0dot/omegatgt)*sin(omegatgt*tf)-(3*x0+2*y0dot/omegatgt)...
    *cos(omegatgt*tf)+(4*x0+2*y0dot/omegatgt);
y=(6*x0+4*y0dot/omegatgt)*sin(omegatgt*tf)+(2*x0dot)/omegatgt...
    *cos(omegatgt*tf)-(6*omegatgt*x0+3*y0dot)*tf+(y0-2*x0dot/omegatgt);
z=z0*cos(omegatgt*tf)+z0dot/omegatgt*sin(omegatgt*tf);

xdot=x0dot*cos(omegatgt*tf)+(3*omegatgt*x0+2*y0dot)*sin(omegatgt*tf);
ydot=(6*omegatgt*x0+4*y0dot)*cos(omegatgt*tf)-2*x0dot*sin(omegatgt*tf)...
    -(6*omegatgt*x0+3*y0dot);
zdot=-z0*omegatgt*sin(omegatgt*tf)+z0dot*cos(omegatgt*tf);

tf=tf/60;

% x0=x0*6378.137;
% y0=y0*6378.137;
% z0=z0*6378.137;
% rtgt=rtgt*6378.137;
% x0dot=x0dot*7.9053838;
% y0dot=y0dot*7.9053838;
% z0dot=z0dot*7.9053838;
% tf=tf*13.446849;
end

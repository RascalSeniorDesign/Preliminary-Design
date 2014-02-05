function [x0dot,y0dot,z0dot] = CWDocking(x0,y0,z0,rtgt,tf)
%The CWDocking equation takes in initial position values in the relative
%target frame, a target oribt, and a target time, and calculates the
%initial relative velocity necessary to acheive docking with a target
%spacecraft. This formulation applies for nearly circular orbits.
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
%Initial Release, CWDocking.m, Tom Moline, 2/5/2014


%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================

omegatgt=sqrt(398600.5/rtgt^3);

y0dot=((6*x0*(omegatgt*tf-sin(omegatgt*tf))-y0)*omegatgt*sind(omegatgt*tf...
    )-2*omegatgt*x0*(4-3*cos(omegatgt*tf))*(1-cos(omegatgt*tf)))...
    /((4*sin(omegatgt*tf)-3*omegatgt*tf)*sin(omegatgt*tf)+4*...
(1-cos(omegatgt*tf))^2);
x0dot=-(omegatgt*x0*(4-3*cos(omegatgt*tf))+2*(1-cos(omegatgt*tf))*y0dot)...
    /sin(omegatgt*tf);
z0dot=-z0*omegatgt*cot(omegatgt*tf);

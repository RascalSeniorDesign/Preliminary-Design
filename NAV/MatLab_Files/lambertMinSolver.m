function [tm,tm_days,am,betam] = lambertMinSolver(r0_,r_)
%The lambertMinSolver function takes in starting and ending position 
%vectors(r0_ and r_). It outputs the minimum time, semi-major axis, and
%eccentic anamoly necessary to make such a transfer
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      r0_       Starting positon vector      3-vector            km
%      r_        Ending position vector       3-vector            km
%      tm        Time of transfer              Scalar       Canonical Time
%      tm_days   Time of Transfer              Scalar            Days
%      am        Min. semi-major axis          Scalar             km
%      betam     Min. Eccentric Anamoly        Scalar            rads
%==========================================================================
%Initial Release, lamberSolver.m, Tom Moline

r0=sqrt(r0_(1)^2+r0_(2)^2+r0_(3)^2); %Positon vector magnitudes, km
r=sqrt(r_(1)^2+r_(2)^2+r_(3)^2);
mu=3986000; %Earth gravitational parameter, km^2/s^3
rearth=6378.1; %Earth Radius, km

c=sqrt(r0^2+r^2-2*r0*r*cosd(theta)); %Chord length, km
s=(r+r0+c)/2; %Semi-perimeter, km
am=s/2; %Minimum semi-major axis, km
betam=2*asin(sqrt((s-c)/s)); %rad
tm=(s^3/8)^(1/2)*(pi-betam+sin(betam)); %canonical time units
TU=sqrt((r0*rearth)^3/mu);
tm_days=tm*TU*(1/(60*60*24)); %Orbit time in days
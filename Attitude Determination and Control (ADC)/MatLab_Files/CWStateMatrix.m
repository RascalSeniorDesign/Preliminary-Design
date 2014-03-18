function [M,N,S,T] = CWStateMatrix(rtgt0_,t)
%The CWStateMatrix function takes in the inertial position of a target
%spacecraft and a time and outputs the state transition matrix that defines
%the relative position and velocity between a target and chaser spacecraft
%at said time.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      rtgt0_    Initial Inertial Target Pos  3x1 Vector          km
%      t         Rendezvous Transfer Time      Scalar             s
%==========================================================================
%Initial Release, CWStateMatrix, Tom Moline, 3/16/2014

%Begin Code

%% Initialize Variables

% Convert Inputs to Canonical Units
t=t/806.811;
rtgt0_=rtgt0_./6378.137;
% Find target orbit magnitude, n, s, and c
rtgt0=sqrt(sum(abs(rtgt0_).^2));
n=sqrt(1/rtgt0^3);
s=sin(n*t);
c=cos(n*t);

%% Create State Transition Matrix Elements

% The State Transition Matrix is Defined as: Phi = [M N; S T]
 M=[ 4-3.*c,     0,       0;...
     6*(s-n*t),  1,       0;...
     0,          0,       c];
 
 S=[3*n*s,        0,     0;...
    -6*n*(1-c),   0,     0;...
    0,            0,     -n*s];

 N=[  s/n,          (2/n)*(1-c),       0;...
     (-2/n)*(1-c),  (4*s-3*n*t)/n,     0;...
      0,             0,                 s/n];
  
 T=[c,      2*s,        0;...
    -2*s,   4*c-3,      0;...
    0,      0,          c];
function [dr_,dv_] = CWPrussing(dr0_,dv0_,rtgt0_,t)
%The CWPrussing function takes in the initial displacement and velocity
%vectors of two spacecraft, as well as the inertial position of a target
%spacecraft and time of transfer, and finds the deltaV necessary for
%rendezvous.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      dr0_      Initial Relative Positoin    3x1 Vector          km
%      dv0_      Initial Relative Velocity    3x1 Vector         km/s
%      rtgt0_    Initial Inertial Target Pos  3x1 Vector          km
%      t         Rendezvous Transfer Time      Scalar             s
%==========================================================================
%Initial Release, CWPrussing, Tom Moline, 3/2/2014


%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
t=t/806.811;
rtgt0_=rtgt0_./6378.137;
dr0_=dr0_./6378.137;
dv0_=dv0_./7.9053838;
rtgt0=sqrt(sum(abs(rtgt0_).^2));
n=sqrt(1/rtgt0^3);
s=sin(n*t);
c=cos(n*t);

%==========================================================================
%                       Define Transition Matrix
%==========================================================================

phi=[4-3*c 0 0 s/n (2/n)*(1-c) 0;...
     6*(s-n*t) 1 0 (-2/n)*(1-c) (4*s-3*n*t)/n 0;...
     0 0 c 0 0 s/n;...
     3*n*s 0 0 c 2*s 0;
     -6*n*(1-c) 0 0 -2*s 4*c-3 0; 0 0 -n*s 0 0 c];
 
 M=phi(1:3,1:3);
 N=phi(1:3, 4:6);
 S=phi(4:6, 1:3);
 T=phi(4:6, 4:6);
 
%==========================================================================
%                       Find Delta dr, dv
%==========================================================================

dr_=M*dr0_+N*dv0_;
dv_=S*dr0_+T*dv0_;
dr_=dr_*6378.137;
dv_=dv_*7.9053838;

% deltaVi=(-inv(N)*M)*dr0_-dv0_;
% deltaVf=((T\N)*M-S)*dr0_;
 
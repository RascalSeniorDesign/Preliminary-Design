function [rp,collision] = hitEarthChecker(rint_,rtgt_,vtransa_,vtransb_,a)
%The hitEarthChecker function takes in the initial and target position
%vectors of a transfer orbit(rint_, rtgt_) and, along with the velocity at
%each point and the semi-major axis of the transfer path, checks if said
%path will result in a collision with the Earth. A +1 indicates no
%collision, a -1 indicates a collision.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      r0_       Starting positon vector      3-vector            km
%      v0_       Starting Velocity vector     3-vector           km/s
%      tf        Flight Time                   Scalar             min
%      r_        Ending Position Vector       3-vector            km
%      v_        Ending Velocity Vector       3-vector           km/s
%==========================================================================
%Initial Release, keplarSolver.m, Tom Moline, 1/31/2014

%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================
rint_=rint_./6378.1;
rtgt_=rtgt_./6378.1;
vtransa_=vtransa_./7.9053838;
vtransb_=vtransb_./7.9053838;
a=a/6378.1;

if dot(rint_,vtransa_)<0.0 && dot(rtgt_,vtransb_)>0.0
    zi=-1/(2*a);
    ht_=cross(rint_,vtransa_);
    ht=sqrt(sum(abs(ht_)).^2);
    p=ht^2;
    e=sqrt((a-p)/a);
    rp=a*(1-e);
    if rp<1
        collision=-1;
    else
        collision=1;
    end
end
        
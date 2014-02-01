function [deltaVa_,deltaVb_] = targetFinder(rint_,rtgt_,vint_,vtgt_,tf)
%The targetFinder function takes in the current position and velocities of
%a interceptor and target spacecraft and a time for which to complete a
%transfer and outputs the deltaV required to accomplish said maneuver.
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      rint_     Interceptor Position vector  3-vector            km
%      rtgt_     Target Position vector       3-vector            km
%      vint_     Interceptr Velocity vector   3-vector            km/s
%      vtgt_     Target Velocity vector       3-vector            km/s
%      tf        Transfer Time                Scalar              s
%==========================================================================
%Initial Release, targetFinder.m, Tom Moline, 1/31/2014

%Begin Code

%==========================================================================
%            Calculate Distance Target Travels During Transfer
%==========================================================================

[rtgtb_,vtgtb_]=keplarSolver(rtgt_,vtgt_,tf);

%==========================================================================
%        Calculate the Transfer Orbit Initial and Final Velocities
%==========================================================================

tm=1; %Taking the Short Path

[vtransa_,vtransb_]=lambertSolver(rint_,rtgtb_,tf,tm);

%==========================================================================
%                          Calculate Delta V
%==========================================================================

deltaVa_=vtransa_-vint_;
deltaVb_=vtgtb_-vtransb_;

%==========================================================================
%                     Find Keps for Transfer Orbit
%==========================================================================

% [p,a,e,i,Omega,omega,theta] = keplarElements(rint_,vtransa_);

%==========================================================================
%           Check if Transfer Orbit Pusts Interceptor Through Earth
%==========================================================================

% [rp,collision]=hitEarthChecker(rint_,rtgt_,vtransa_,vtransb_,a);
% 
% if collision==1
%     disp('You did not hit Earth, congrats!')
% else
%     disp('You crashed an burned')
% end



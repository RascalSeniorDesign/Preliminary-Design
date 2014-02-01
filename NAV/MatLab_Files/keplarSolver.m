function [r_,v_] = keplarSolver(r0_,v0_,tf)
%The keplarSolver function takes in position and velcity vectors at a given
%time (r0_, v0_) and, provided an end time (tf), finds the position and
%velocity vectors at some other point on the spacecraft orbit or
%trajectory. This works for all types of orbits/trajectories and makes use
%of universal variables and canonical time in order to speed up operations.
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
r0_=r0_./6378.1; %Convert from km to canonical distance, DU
v0_=v0_./7.9053838; %Convert from km/s to DU/TU
r0=sqrt(sum(abs(r0_).^2)); %Initial displacement magnitude
v0=sqrt(sum(abs(v0_).^2)); %Initial velocity magnitude
tf=tf/13.446849;
%==========================================================================
%                       Define Universal Variables
%==========================================================================
zi=v0^2/2-1/r0;
a=-1/(2*zi);
alpha=1/a; %Universal variable, replaces eccentricity, DU^-1

if alpha>0.000001 %For circular or elliptical case
    xi0=tf*alpha;
elseif alpha==1 %For case when guess is too close to solution to converge
    xi0=tf*alpha*.95;
elseif abs(alpha)<0.000001 %For parabolic case
    h_=cross(r0_,v0_);
    h=sqrt(sum(abs(h_).^2));
    p=h^2;
    s=acotd(3*sqrt(1/p^3)*tf*.5);
    w=atand((tand(s))^1/3);
    xi0=sqrt(p)*2*cotd(2*w);
elseif alpha<-0.000001 %For hyperbolic case
    a=1/alpha;
    xi0=sign(tf)*sqrt(-a)*log((-2*alpha*tf)/(dot(r0_,v0_)+sign(tf)*sqrt(-a)...
        *(1-r0*alpha)));
end

%==========================================================================
%                      Set Up and Run Iteration Loop
%==========================================================================
i=1;
Nmax=10^3;
dotrv0=dot(r0_,v0_); %Prevents doing dot product for each iteration

while i<Nmax
    psi=xi0^2*alpha;
    if psi>10^-6 %Checks if elliptical
        c2=(1-cos(sqrt(psi)))/psi;
        c3=(sqrt(psi)-sin(sqrt(psi)))/sqrt(psi^3);
    elseif psi<-10^-6 %Checks if hyperbolic
        c2=(1-cosh(sqrt(-psi)))/psi;
        c3=(sqrt(-psi)-sinh(sqrt(-psi)))/sqrt((-psi)^3);
    else %If equal to zero
        c2=.5;
        c3=1/6;
    end
    r=xi0^2*c2+dotrv0*xi0*(1-psi*c3)+r0*(1-psi*c2);
    xi1=xi0+(tf-xi0^3*c3-dotrv0*xi0^2*c2-r0*xi0*(1-psi*c3))/r;
    if abs(xi1-xi0)<10^-6
        break
    else
        xi0=xi1;
        i=i+1;
    end
end

%==========================================================================
%                               Find r_ and v_
%==========================================================================
f=1-(xi1^2/r0)*c2;
g=tf-xi1^3*c3;
gdot=1-(xi1^2/r)*c2;
fdot=(1/(r*r0))*xi1*(psi*c3-1);
r_=f.*r0_+g.*v0_;
v_=fdot.*r0_+gdot.*v0_;
r_=r_.*6378.1; %Converting final position from canonical to km
v_=v_.*7.9053838; %Converting final velocity from canonical to km/s
    
        
    
    
    








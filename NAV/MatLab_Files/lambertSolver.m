function [v0_,v_] = lambertSolver(r0_,r_,tf,tm)
%The lambertSolver function takes in starting and ending position vectors
%(r0_ and r_), as well as the time to reach an orbit that intersects each
%of those two points (tf), and the manner in which one wishes to do
%so(Long, +1, or short, -1, as specified by tm). It outputs the intial and
%final velocity vectors necessary to reach said orbit through the use of
%Universal Variable Algorithm, as described in Vallado
%
%==========================================================================
% Variable Name  Variable Description      Variable Type    Variable Units
%==========================================================================
%      r0_       Starting positon vector      3-vector            km
%      r_        Ending position vector       3-vector            km
%      tf        Time of transfer              Scalar       Canonical Time
%      tm        Method of transfer            Scalar           1 or -1
%      v0_       Starting velocity vector     3-vector           km/s
%      v_        Ending velocity vector       3-vector           km/s
%==========================================================================
%Initial Release, lamberSolver.m, Tom Moline

%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================

r0_=r0_./6378.1;
r_=r_./6378.1;

r0=sqrt(sum(abs(r0_).^2)); %Positon vector magnitudes, km
r=sqrt(sum(abs(r_).^2));
mu=398600; %Earth gravitational parameter, km^2/s^3
theta=acosd(dot((r0_),(r_))/(r0*r)); %True Anomaly, deg
tf=tf/13.446849;
%==========================================================================
%                    Set Up Universal Variables
%============================== ============================================
A=tm*sqrt(r*r0*(1+cosd(theta)));
c2=.5;
c3=1/6;
phin=0.0;
phiup=4*pi^2;
philow=0;
i=1;
Nmax=10^3;
%==========================================================================
%                    Run Univeral Variable Algorithm
%==========================================================================
while i<Nmax
    yn=r+r0+(A*(phin*c3-1)/sqrt(c2)); 
    if A>0.0 && yn<0.0
        while yn<0.0
            philow=philow+.01;
            yn=r0+r+(A*(philow*c3-1)/sqrt(c2));
        end
    end
    xin=sqrt(yn/c2);
    tn=(xin^3*c3+A*sqrt(yn))/sqrt(1);
    if tn<tf
        philow=phin;
    else
        phiup=phin;
    end
    phintemp=(phiup+philow)/2;
    if phintemp>10^-6
        c2=(1-cos(sqrt(phintemp)))/phintemp;
        c3=(sqrt(phintemp)-sin(sqrt(phintemp)))/sqrt(phintemp^3);
    elseif phintemp<-10^-6
        c2=(1-cosh(sqrt(-phintemp)))/phintemp;
        c3=(sqrt(-phintemp)-sinh(sqrt(-phintemp)))/sqrt((-phintemp)^3);
    else
        c2=.5;
        c3=1/6;
    end
    phin=phintemp;
    if abs(tn-tf)<10^-6
        break
    else
        i=i+1;
    end
end
%==========================================================================
%                     Find Velocity Vectors
%==========================================================================
f=1-yn/r0;
gdot=1-yn/r;
g=A*sqrt(yn/1);
v0_=(r_-f*r0_)./g;
v_=(gdot*r_-r0_)./g;

v0_=v0_.*7.9053838;
v_=v_.*7.9053838;
end









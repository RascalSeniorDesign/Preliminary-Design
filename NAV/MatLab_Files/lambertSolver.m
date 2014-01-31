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
%      tf        Time of transfer              Scalar             min
%      tm        Method of transfer            Scalar           1 or -1
%      v0_       Starting velocity vector     3-vector           km/s
%      v_        Ending velocity vector       3-vector           km/s
%==========================================================================
%Initial Release, lamberSolver.m, Tom Moline, 1/29/2014
%Updated and functional, lambertSolver.m, Tom Moline, 1/31/2014

%Begin Code

%==========================================================================
%                       Initialize Variables
%==========================================================================

r0_=r0_./6378.1; %Convert from km to Distance Units (related to earth radi)
r_=r_./6378.1;

r0=sqrt(sum(abs(r0_).^2)); %Positon vector magnitudes, DU
r=sqrt(sum(abs(r_).^2));
theta=acosd(dot((r0_),(r_))/(r0*r)); %True Anomaly, deg
tf=tf/13.446849; %Canonical time, TU
%==========================================================================
%                    Set Up Universal Variables
%============================== ============================================
A=tm*sqrt(r*r0*(1+cosd(theta)));
c2=.5;
c3=1/6;
phin=0.0;
phiup=4*pi^2; %Related to square of eccentric anamoly, rad^2
philow=0;
i=1; %Sitting up iterative while loop
Nmax=10^3; %Max number of iterations
%==========================================================================
%     Run Univeral Variable Algorithm (Algorithm 57, Vallado, Page 488)
%==========================================================================
while i<Nmax
    yn=r+r0+(A*(phin*c3-1)/sqrt(c2));%Relative area realted to current phin
    if A>0.0 && yn<0.0
        while yn<0.0
            philow=philow+.01;
            yn=r0+r+(A*(philow*c3-1)/sqrt(c2));
        end
    end
    xin=sqrt(yn/c2);
    tn=(xin^3*c3+A*sqrt(yn));
    if tn<tf %Bisection method
        philow=phin;
    else
        phiup=phin;
    end
    phintemp=(phiup+philow)/2;
    if phintemp>10^-6 %Checks if elliptical
        c2=(1-cos(sqrt(phintemp)))/phintemp;
        c3=(sqrt(phintemp)-sin(sqrt(phintemp)))/sqrt(phintemp^3);
    elseif phintemp<-10^-6 %Checks if hyperbolic
        c2=(1-cosh(sqrt(-phintemp)))/phintemp;
        c3=(sqrt(-phintemp)-sinh(sqrt(-phintemp)))/sqrt((-phintemp)^3);
    else %If equal to zero
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
g=A*sqrt(yn);
v0_=(r_-f*r0_)./g; %Velocity at initial point, DU/TU
v_=(gdot*r_-r0_)./g; %Velocity at final point, DU/TU

v0_=v0_.*7.9053838; %Velocity at initial point, km/s
v_=v_.*7.9053838; %Velocity at final point, km/s
end









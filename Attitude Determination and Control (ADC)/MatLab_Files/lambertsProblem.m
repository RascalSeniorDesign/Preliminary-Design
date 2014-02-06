r1=1; %AU
r2=1.524; %AU
theta=75; %degrees
mu=398600;

c=sqrt(r1^2+r2^2-2*r1*r2*cosd(theta));
s=(r1+r2+c)/2;

am=s/2;
betam=2*asin(sqrt((s-c)/s));
tm=(1/sqrt(1))*(s^3/8)^(1/3)*(pi-betam+sin(betam));


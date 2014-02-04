tf=linspace(0,250,1000);

x=zeros(1,length(tf));
y=zeros(1,length(tf));
z=zeros(1,length(tf));
xdot=zeros(1,length(tf));
ydot=zeros(1,length(tf));
zdot=zeros(1,length(tf));

rtgt=6378.137+580;
x0=0;
y0=0;
z0=0;
x0dot=-.1;
y0dot=-.04;
z0dot=-.02;

for i=1:length(tf)
    [x(i),y(i),z(i),xdot(i),ydot(i),zdot(i)] = CWSolver(x0,y0,z0,x0dot,y0dot,z0dot,rtgt,tf(i));
end

figure(1)
plot(tf,x,tf,y,tf,z)
xlabel('Time (min)')
ylabel('Relative Position (m)')
legend('x','y','z')

figure(2)
plot(tf,xdot,tf,ydot,tf,zdot)
xlabel('Time (min)')
ylabel('Relative Velocity (m/s)')
legend('dx','dy','dz')
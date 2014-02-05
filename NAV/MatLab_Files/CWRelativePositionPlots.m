clear
clc

tf=linspace(5,250,100);

x=zeros(length(tf),length(tf));
y=zeros(length(tf),length(tf));
z=zeros(length(tf),length(tf));
xdot=zeros(length(tf),length(tf));
ydot=zeros(length(tf),length(tf));
zdot=zeros(length(tf),length(tf));
deltaV=zeros(length(tf),length(tf));

rtgt=6378.137+580;
x0=-100;
y0=-100;
z0=-100;
% x0dot=-.1;
% y0dot=-.04;
% z0dot=-.02;

for i=1:length(tf)
    [x0dot,y0dot,z0dot] = CWDocking(x0,y0,z0,rtgt,tf(i));
    for j=1:length(tf)
        [x(i,j),y(i,j),z(i,j),xdot(i,j),ydot(i,j),zdot(i,j)] = CWSolver(x0,y0,z0,x0dot,y0dot,z0dot,rtgt,tf(j));
    end
end

% figure (1)
% plot(tf,deltaV,tf,abs(xdot),tf,abs(ydot),tf,abs(zdot))
% xlabel('Time of Transfer (min)')
% ylabel('Delta V Required (min)')

% figure(2)
% plot(tf,x,tf,y,tf,z)
% xlabel('Time (min)')
% ylabel('Relative Position (m)')
% legend('x','y','z')
% 
% figure(3)
% plot(tf,xdot,tf,ydot,tf,zdot)
% xlabel('Time (min)')
% ylabel('Relative Velocity (m/s)')
% legend('dx','dy','dz')

figure(1)
mesh(tf,tf,y,'EdgeColor','black','FaceColor','None')
grid on
set(gca,'GridLineStyle','-')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
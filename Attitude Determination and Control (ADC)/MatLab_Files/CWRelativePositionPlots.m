clear
clc

tf=linspace(0,250,100);
tp=linspace(5,250,100);

x=zeros(length(tf),length(tp));
y=zeros(length(tf),length(tp));
z=zeros(length(tf),length(tp));
xdot=zeros(length(tf),length(tp));
ydot=zeros(length(tf),length(tp));
zdot=zeros(length(tf),length(tp));
deltaV=zeros(length(tf),length(tp));
x0dot=zeros(1,length(tp));
y0dot=zeros(1,length(tp));
z0dot=zeros(1,length(tp));


rtgt=6378.137+580;
x0=1000000;
y0=1000000;
z0=1000000;
% x0dot=-.1;
% y0dot=-.04;
% z0dot=-.02;

for i=1:length(tf)
    [x0dot(i),y0dot(i),z0dot(i)] = CWDocking(x0,y0,z0,rtgt,tf(i));
    deltaV(i)=sqrt(abs(x0dot(i))^2+abs(y0dot(i))^2+abs(z0dot(i))^2);
    for j=1:length(tp)
        [x(i,j),y(i,j),z(i,j),xdot(i,j),ydot(i,j),zdot(i,j)] = CWSolver(x0,y0,z0,x0dot(i),y0dot(i),z0dot(i),rtgt,tp(j));
    end
end

figure (1)
plot(tf,deltaV,tf,x0dot,tf,y0dot,tf,z0dot)
xlabel('Time of Transfer (min)')
ylabel('Delta V Required (min)')

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

figure(2)
subplot(2,1,1)
mesh(tf,tp,x,'EdgeColor','black','FaceColor','None')
view(44,28)
grid on
set(gca,'GridLineStyle','-')
subplot(2,1,2)
mesh(tf,tp,y,'EdgeColor','black','FaceColor','None')
view(44,28)
grid on
set(gca,'GridLineStyle','-')

figure(3)
mesh(tf,tp,z,'EdgeColor','black','FaceColor','None')
view(44,28)
grid on
set(gca,'GridLineStyle','-')
xlabel('Transfer Time (m)')
ylabel('Plot Time (m)')
zlabel('Relative Y Position (m)')
clear;
clc;

function rot=R(fi)
    rot=[cos(fi), -sin(fi); sin(fi), cos(fi)];
endfunction

T=10;
h=0.01;
t=[0:h:T];
n=length(t)

#Stałe:
global I=[1, 0; 0, 1];
global Q=[0, -1; 1, 0]; 

global l1=1;
global m1=1;
global J1=1/12 * m1 * l1^2;
global l2=1;
global m2=1;
global J2=1/12 * m2 * l2^2;
global l3=1;
global m3=1;
global J3=1/12 * m3 * l3^2;

global M=diag([m1;m1;J1;m2;m2;J2;m3;m3;J3]);

global s_A1=[-l1/2;0];
global s_B1=[l1/2;0];
global s_B2=[-l2/2;0];
global s_C2=[l2/2;0];
global s_C3=[-l3/2;0];
global s_D3=[l3/2;0];

l=[l1,l2,l3];

global S1=zeros(2,3);
global S2=zeros(2,3);
for part=[1:1:3]
      S1(:,part)=[-l(part)/2;0];
      S2(:,part)=[l(part)/2;0];
endfor


global grav=-9.81;

theta1=0;
theta2=0;
theta3=0;
thetadot1=0;
thetadot2=0;
thetadot3=0;

x0=[theta1,theta2, theta3, thetadot1,thetadot2, thetadot3];


[x, istate, msg] = lsode (@three_parts_dynamics_function, x0, t);


fi1=x(:,1);
fi2=x(:,2)+x(:,1);
fi3=x(:,3)+x(:,2)+x(:,1);
fidot1=x(:,4);
fidot2=x(:,4)+x(:,5);
fidot3=x(:,4)+x(:,5)+x(:,6);
x3=l1*cos(fi1)+l2*cos(fi2)+l3/2*cos(fi3);
y3=l1*sin(fi1)+l2*sin(fi2)+l3/2*sin(fi3);

load t.txt
load x3.txt
load y3.txt

#plot(t(1:1:401),x3(1:1:401),"LineWidth",3,t(1:1:401),y3(1:1:401),"LineWidth",3,t_free(1:1:401),x_free(1:1:401),"LineWidth",2,t_free(1:1:401),y_free(1:1:401),"LineWidth",2)
#legend("3.part x - Octave", "3.part y - Octave","3.part x - FreeDyn","3.part y - FreeDyn")
#xlabel("czas[s]")
#ylabel("położenie[m]")
#set(gca, "fontsize", 20)
#plot(xc,yc,xks,yks)
#plot([zeros(1,101);xb],[zeros(1,101);yb],[xb;xc],[yb;yc])

plot(t(1:1:401),abs(x3(1:1:401)-x_free(1:1:401)),"LineWidth",3,t(1:1:401),abs(y3(1:1:401)-y_free(1:1:401)),"LineWidth",3)
set(gca, "fontsize", 20)
legend("moduł różnicy x", "moduł różnicy y","location","northwest")
xlabel("czas[s]")

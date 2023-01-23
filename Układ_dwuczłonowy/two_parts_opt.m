clear;
clc;

function rot=R(fi)
    rot=[cos(fi), -sin(fi); sin(fi), cos(fi)];
endfunction

T=1;
global h=0.05;
t=[0:h:T];
global n=length(t)

#Stałe:
global I=[1, 0; 0, 1];
global Q=[0, -1; 1, 0]; 

global l1=1;
global m1=1;
global J1=1;
#global J1=1/12 * m1 * l1^2;
global l2=1;
global m2=1;
global J2=1;
#global J2=1/12 * m2 * l2^2;

global s_A1=[-l1/2;0];
global s_B1=[l1/2;0];
global s_B2=[-l2/2;0];

global M=diag([m1;m1;J1;m2;m2;J2]);

l=[l1,l2];

global S1=zeros(2,2);
global S2=zeros(2,2);
for part=[1:1:2]
      S1(:,part)=[-l(part)/2;0];
      S2(:,part)=[l(part)/2;0];
endfor

global grav=0;
#g=[0;grav*m1;0;0;grav*m2;0];
global xp=2
global yp=0
global xk=0.5
global yk=1.5

#theta1=0;
#theta2=0;
#thetadot1=0;
#thetadot2=0;

load x.mat
k = rows(x)/6;

tp=[0:T/(k-1):T]';

p10 = interp1(tp, x([1:k]), t');
p20 = interp1(tp, x([k+1:2*k]), t');
v10 = interp1(tp, x([2*k+1:3*k]), t');
v20 = interp1(tp, x([3*k+1:4*k]), t');
u10 = interp1(tp, x([4*k+1:5*k]), t');
u20 = interp1(tp, x([5*k+1:6*k]), t');
clear x
#x0 = [p10; p20; v10; v20; u10; u20];
x0=ones(6*n,1);
#p10 = t'/T *2.11;
#p20 = t'/T *(-1.55);
#p10 = 0.1+t'/T *0.6;
#p20 = -0.3-t'/T *0.2;
#x0 = [p10; p20; ones(2*n,1); zeros(2*n,1)];
#x0 = [p10; p20; zeros(4*n,1)];

function r = g(x)
    global n;
    global h;
    global m1;
    global m2;
    global l1;
    global l2;
    global grav;
    global xk;
    global yk;
    global xp;
    global yp;
    global S_1
    global S_2

    theta1 = x([1:n]);
    theta2 = x([n+1:2*n]);
    thetadot1 = x([2*n+1:3*n]);
    thetadot2 = x([3*n+1:4*n]);
    q=[theta1';theta2';thetadot1';thetadot2'];
    u=[x([4*n+1:5*n])';x([5*n+1:6*n])'];

    #r =[thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n); l1*cos(theta1(n))+l2*cos(theta2(n)+theta1(n))-xk; l1*sin(theta1(n))+l2*sin(theta2(n)+theta1(n))-yk; l1*cos(theta1(1))+l2*cos(theta2(1)+theta1(1))-xp; l1*sin(theta1(1))+l2*sin(theta2(1)+theta1(1))-yp]; # Warunki brzegowe
    r =[theta1(1); theta1(n)-pi/2; theta2(1); theta2(n)-pi/2; thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n)];
    for i=[1:n-1]
        vec = q(:,i+1)-q(:,i)-1/2*h*(two_parts_dynamics_function(q(:,i+1),u(:,i+1))+two_parts_dynamics_function(q(:,i),u(:,i)));
        r=[r;vec];
    endfor

endfunction

function obj = phi(x)
    global n;
    global h;
    obj = 0;
    u1 = x([4*n+1:5*n]);
    u2 = x([5*n+1:6*n]);

    for i=[1:n-1]
        obj = obj + 1/2 * h * (u1(i)^2 + u1(i+1)^2 + u2(i)^2 + u2(i+1)^2);
    endfor

endfunction
#sqp (x0, @phi, @g, [],[],[],300)
tic();
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, [],[],[],300);
time_of_execution=toc()

obj
info
iter
time_of_execution

theta1 = x([1:n]);
theta2 = x([n+1:2*n]);
thetadot1 = x([2*n+1:3*n]);
thetadot2 = x([3*n+1:4*n]);
u1 = x([4*n+1:5*n]);
u2 = x([5*n+1:6*n]);
#save wyniki.mat theta1 theta2 thetadot1 thetadot2 u1 u2
#save x.mat x
#plot(t,theta1,t,theta2);
#plot(t,u1,t,u2)
save two_parts_por.mat

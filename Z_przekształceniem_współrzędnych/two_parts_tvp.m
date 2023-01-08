clear;
clc;

function rot=R(fi)
    rot=[cos(fi), -sin(fi); sin(fi), cos(fi)];
endfunction

T=1;
global h=0.1;
t=[0:h:T];
global n=length(t)

#Stałe:
global I=[1, 0; 0, 1];
global Q=[0, -1; 1, 0]; 

global l1=1;
global m1=1;
global J1=1/12 * m1 * l1^2;
global l2=1.2;
global m2=1;
global J2=1/12 * m2 * l2^2;

global M=diag([m1;m1;J1;m2;m2;J2]);

l=[l1,l2];

global S1=zeros(2,2);
global S2=zeros(2,2);
for part=[1:1:2]
      S1(:,part)=[-l(part)/2;0];
      S2(:,part)=[l(part)/2;0];
endfor

#global theta1k=0.7;
#global theta2k=-0.5;
#global theta1p=0.1;
#global theta2p=-0.3;
global xp=2;
global yp=0;
global xk=0.5;
global yk=1.5;

global grav=-9.81;
#g=[0;grav*m1;0;0;grav*m2;0];

[xkon1,xkondot1,xkondotdot1]=function_tvp(T,n,xp,xk);
[ykon1,ykondot1,ykondotdot1]=function_tvp(T,n,yp,yk);
global xkon=xkon1;
global ykon=ykon1;
#[theta2,thetadot2,thetadotdot2]=function_tvp(T,n,theta2p,theta2k);
#x0=[theta1;theta2;thetadot1;thetadot2;thetadotdot1;thetadotdot2];
x0=ones(6*n,1);

function r = g(x)
    global n;
    global h;
    global m1;
    global m2;
    global l1;
    global l2;
    global grav;
    global theta1k;
    global theta2k;
    global theta1p;
    global theta2p;
    global xp;
    global yp;
    global xk;
    global yk;

    theta1 = x([1:n]);
    theta2 = x([n+1:2*n]);
    thetadot1 = x([2*n+1:3*n]);
    thetadot2 = x([3*n+1:4*n]);
    q=[theta1';theta2';thetadot1';thetadot2'];
    u=[x([4*n+1:5*n])';x([5*n+1:6*n])'];

    r =[thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n); l1*cos(theta1(n))+l2*cos(theta2(n)+theta1(n))-xk; l1*sin(theta1(n))+l2*sin(theta2(n)+theta1(n))-yk; l1*cos(theta1(1))+l2*cos(theta2(1)+theta1(1))-xp; l1*sin(theta1(1))+l2*sin(theta2(1)+theta1(1))-yp]; # Warunki brzegowe
    #r =[theta1(1)-theta1p; theta1(n)-theta1k; theta2(1)-theta2p; theta2(n)-theta2k; thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n)];
    for i=[1:n-1]
        vec = q(:,i+1)-q(:,i)-1/2*h*(two_parts_dynamics_function(q(:,i+1),u(:,i+1))+two_parts_dynamics_function(q(:,i),u(:,i)));
        r=[r;vec];
    endfor

endfunction

function obj = phi(x)
    global n;
    global h;
    global l1;
    global l2;
    global xkon;
    global ykon;
    obj = 0;
    #u1 = x([4*n+1:5*n]);
    #u2 = x([5*n+1:6*n]);

    #for i=[1:n-1]
        #obj = obj + 1/2 * h * (u1(i)^2 + u1(i+1)^2 + u2(i)^2 + u2(i+1)^2);
    #endfor
    theta1 = x([1:n]);
    theta2 = x([n+1:2*n]);
    for i=[1:n-1]
        obj = obj + (l1*cos(theta1(i))+l2*cos(theta2(i)+theta1(i))-xkon(i))^2+(l1*sin(theta1(i))+l2*sin(theta2(i)+theta1(i))-ykon(i))^2;
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
save wyniki.mat theta1 theta2 thetadot1 thetadot2 u1 u2
save x.mat x
#plot(t,theta1,t,theta2+theta1);
#plot(t,u1,t,u2)
xb=l1*cos(theta1);
yb=l1*sin(theta1);
xc=l1*cos(theta1)+l2*cos(theta2+theta1);
yc=l1*sin(theta1)+l2*sin(theta2+theta1);
plot(l1*cos(theta1)+l2*cos(theta2+theta1), l1*sin(theta1)+l2*sin(theta2+theta1), "LineWidth",3)
hold
plot([zeros(1,n);xb'],[zeros(1,n);yb'],[xb';xc'],[yb';yc'])
axis("equal")
clear;
clc;

tic();

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
global l2=1.1;
global m2=1;
global J2=1/12 * m2 * l2^2;
global l3=1.2;
global m3=1;
global J3=1/12 * m3 * l3^2;

global M=diag([m1;m1;J1;m2;m2;J2;m3;m3;J3]);
global s_A1=[-l1/2;0];
global s_B1=[l1/2;0];
global s_B2=[-l2/2;0];
global s_C2=[l2/2;0];
global s_C3=[-l3/2;0];
global s_D3=[l3/2;0];
global grav=-9.81;
global xp=2
global yp=0.5
global xk=2.5
global yk=1

l=[l1,l2,l3];

global S1=zeros(2,3);
global S2=zeros(2,3);
for part=[1:1:3]
      S1(:,part)=[-l(part)/2;0];
      S2(:,part)=[l(part)/2;0];
endfor


x0=zeros(9*n,1);
[xkon1,xkondot1,xkondotdot1]=function_tvp(T,n,xp,xk);
[ykon1,ykondot1,ykondotdot1]=function_tvp(T,n,yp,yk);
global xkon=xkon1;
global ykon=ykon1;


function r = g(x)
    global n;
    global h;
    global m1;
    global m2;
    global m3;
    global l1;
    global l2;
    global l3;
    global grav;
    global xk;
    global yk;
    global xp;
    global yp;

    theta1 = x([1:n]);
    theta2 = x([n+1:2*n]);
    theta3 = x([2*n+1:3*n]);
    thetadot1 = x([3*n+1:4*n]);
    thetadot2 = x([4*n+1:5*n]);
    thetadot3 = x([5*n+1:6*n]);
    u = [x([6*n+1:7*n])';x([7*n+1:8*n])';x([8*n+1:9*n])'];
    q = [theta1';theta2'; theta3'; thetadot1'; thetadot2'; thetadot3'];

    r =[thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n); thetadot3(1); thetadot3(n); l1*cos(theta1(n))+l2*cos(theta2(n)+theta1(n))+l3*cos(theta3(n)+theta2(n)+theta1(n))-xk; l1*sin(theta1(n))+l2*sin(theta2(n)+theta1(n))+l3*sin(theta3(n)+theta2(n)+theta1(n))-yk; l1*cos(theta1(1))+l2*cos(theta2(1)+theta1(1))+l3*cos(theta3(1)+theta2(1)+theta1(1))-xp; l1*sin(theta1(1))+l2*sin(theta2(1)+theta1(1))+l3*sin(theta3(1)+theta2(1)+theta1(1))-yp]; # Warunki brzegowe

    #r=[theta1(1)+1; theta1(n)-1; theta2(1)+1; theta2(n)-1; theta3(1)+1; theta3(n)-1;
    #thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n); thetadot3(1); thetadot3(n);];

    for i=[1:n-1]
        vec = q(:,i+1)-q(:,i)-1/2*h*(three_parts_dynamic_function(q(:,i+1),u(:,i+1))+three_parts_dynamic_function(q(:,i),u(:,i)));
        r=[r;vec];
    endfor

endfunction

function obj = phi(x)
    global n;
    global h;
    global xkon;
    global ykon;
    global l1;
    global l2;
    global l3;
    obj = 0;
    
    theta1 = x([1:n]);
    theta2 = x([n+1:2*n]);
    theta3 = x([2*n+1:3*n]);
    for i=[1:n-1]
        obj = obj + (l1*cos(theta1(i))+l2*cos(theta2(i)+theta1(i))+l3*cos(theta3(i)+theta2(i)+theta1(i))-xkon(i))^2+(l1*sin(theta1(i))+l2*sin(theta2(i)+theta1(i))+l3*sin(theta3(n)+theta2(n)+theta1(n))-ykon(i))^2;
    endfor

endfunction
#sqp (x0, @phi, @g, [],[],[],300)
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, [],[],[],300);
time_of_execution=toc()

obj
info
iter

theta1 = x([1:n]);
theta2 = x([n+1:2*n]);
theta3 = x([2*n+1:3*n]);
thetadot1 = x([3*n+1:4*n]);
thetadot2 = x([4*n+1:5*n]);
thetadot3 = x([5*n+1:6*n]);
u1 = x([6*n+1:7*n]); 
u2 = x([7*n+1:8*n]);
u3 = x([8*n+1:9*n]);

#plot(t,theta1,t, theta2, t,theta3)
plot(l1*cos(theta1)+l2*cos(theta2+theta1)+l3*cos(theta3+theta2+theta1), l1*sin(theta1)+l2*sin(theta2+theta1)+l3*sin(theta3+theta2+theta1), "LineWidth",3)
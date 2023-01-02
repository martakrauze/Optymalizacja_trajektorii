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

global l1=1.2;
global m1=1;
global J1=1/12 * m1 * l1^2;
global l2=1.2;
global m2=1;
global J2=1/12 * m2 * l2^2;
global l3=1.2;
global m3=1;
global J3=1/12 * m2 * l2^2;

global M=diag([m1;m1;J1;m2;m2;J2;m3;m3;J3]);
global s_A1=[-l1/2;0];
global s_B1=[l1/2;0];
global s_B2=[-l2/2;0];
global s_C2=[l2/2;0];
global s_C3=[-l3/2;0];
global s_D3=[l3/2;0];
global grav=-9.81;
#g=[0;grav*m1;0;0;grav*m2;0];
global xp=2
global yp=0
global xk=-2
global yk=0

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
x0 = [p10; p20; v10; v20; u10; u20];
#x0=zeros(6*n,1);
#p10 = t'/T *2.11;
#p20 = t'/T *(-1.55);
p10 = t'/T *-pi;
p20 = t'/T *(0);
#x0 = [p10; p20; zeros(4*n,1)];



function r=fcn(x,u)
    global I
    global Q 
    global M
    global s_A1
    global s_B1
    global s_B2
    global grav
    global m1
    global m2

    g=[0;grav*m1;u(1)-u(2);0;grav*m2;u(2)];

    fi1=x(1);
    fi2=x(2)+x(1);
    fidot1=x(3);
    fidot2=x(4)+x(3);  
    thetadot=[x(3);x(4)];

    s_A1_g=R(fi1)*s_A1;
    s_B1_g=R(fi1)*s_B1;
    s_B2_g=R(fi2)*s_B2;

    D=[I, Q * s_A1_g, zeros(2,3);
       I, Q * s_B1_g, -I, -Q * s_B2_g];
    

    d_11=-Q * s_A1_g;
    d_21=Q*(-s_A1_g+s_B1_g-s_B2_g);
    d_22=-Q*s_B2_g;

    B=[d_11, zeros(2,1);
          1, 0;
       d_21, d_22;
          1, 1];
    
    d_11_dot=s_A1_g*fidot1;
    d_21_dot=-(-s_A1_g*fidot1+s_B1_g*fidot1-s_B2_g*fidot2);
    d_22_dot=s_B2_g*fidot2;

    Bdot=[d_11_dot, zeros(2,1);
          0, 0;
       d_21_dot, d_22_dot;
          0, 0];

    M_new=B' * M * B;
    g_new=B' * (g - M * Bdot * thetadot);

    r1=x(3);
    r2=x(4);
    r3=M_new\g_new;
    r=[r1;r2;r3];

endfunction


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

    theta1 = x([1:n]);
    theta2 = x([n+1:2*n]);
    thetadot1 = x([2*n+1:3*n]);
    thetadot2 = x([3*n+1:4*n]);
    u1 = x([4*n+1:5*n]);
    u2 = x([5*n+1:6*n]);

    r =[theta1(1); theta2(1); thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n); l1*cos(theta1(n))+l2*cos(theta2(n)+theta1(n))-xk; l1*sin(theta1(n))+l2*sin(theta2(n)+theta1(n))-yk]; # Warunki brzegowe

    for i=[1:n-1]
        q_i=[theta1(i);theta2(i);thetadot1(i);thetadot2(i)];
        q_ip1=[theta1(i+1);theta2(i+1);thetadot1(i+1);thetadot2(i+1)];
        u_i=[u1(i);u2(i)];
        u_ip1=[u1(i+1);u2(i+1)];
        vec = q_ip1-q_i-1/2*h*(fcn(q_ip1,u_ip1)+fcn(q_i,u_i));
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
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, [],[],[],300);

obj
info
iter

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
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


x0=ones(9*n,1);

function r=fcn(x,u)
    global I
    global Q 
    global M
    global s_A1
    global s_B1
    global s_B2
    global s_C2
    global s_C3
    global s_D3
    global grav
    global m1
    global m2
    global m3

    g=[0;grav*m1;u(1)-u(2);0;grav*m2;u(2)-u(3);0;grav*m3;u(3)];

    fi1=x(1);
    fi2=x(2)+x(1);
    fi3=x(3)+x(2)+x(1);
    fidot1=x(4);
    fidot2=x(5)+x(4);
    fidot3=x(6)+x(5)+x(4);  
    thetadot=[x(4);x(5);x(6)];

    s_A1_g=R(fi1)*s_A1;
    s_B1_g=R(fi1)*s_B1;
    s_B2_g=R(fi2)*s_B2;
    s_C2_g=R(fi2)*s_C2;
    s_C3_g=R(fi3)*s_C3;
    s_D3_g=R(fi3)*s_D3;

    D=[I, Q * s_A1_g, zeros(2,3), zeros(2,3);
       I, Q * s_B1_g, -I, -Q * s_B2_g, zeros(2,3);
       zeros(2,3), I, Q * s_C2_g, -I, -Q * s_C3_g];
    

    d_11=-Q * s_A1_g;
    d_21=Q*(-s_A1_g+s_B1_g-s_B2_g);
    d_22=-Q*s_B2_g;
    d_31=Q*(-s_A1_g+s_B1_g-s_B2_g+s_C2_g-s_C3_g);
    d_32=Q*(-s_B2_g+s_C2_g-s_C3_g);
    d_33=-Q*s_C3_g;

    B=[d_11, zeros(2,1), zeros(2,1);
          1,          0,          0;
       d_21,       d_22, zeros(2,1);
          1,          1,          0;
       d_31,       d_32,       d_33;
          1,          1,         1];
    
    d_11_dot=s_A1_g*fidot1;
    d_21_dot=-(-s_A1_g*fidot1+s_B1_g*fidot1-s_B2_g*fidot2);
    d_22_dot=s_B2_g*fidot2;
    d_31_dot=-(-s_A1_g*fidot1+s_B1_g*fidot1-s_B2_g*fidot2+s_C2_g*fidot2-s_C3_g*fidot3);
    d_32_dot=-(-s_B2_g*fidot2+s_C2_g*fidot2-s_C3_g*fidot3);
    d_33_dot=s_C3_g*fidot3;

    Bdot=[d_11_dot, zeros(2,1), zeros(2,1);
                 0,          0,          0;
          d_21_dot,   d_22_dot, zeros(2,1);
                 0,          0,          0;
          d_31_dot,   d_32_dot,   d_33_dot,
                 0,          0,          0];
        

    M_new=B' * M * B;
    g_new=B' * (g - M * Bdot * thetadot);

    r1=x(4);
    r2=x(5);
    r3=x(6);
    r4=M_new\g_new;
    r=[r1;r2;r3;r4];

endfunction


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
    obj = 0;
    u1 = x([6*n+1:7*n]);
    u2 = x([7*n+1:8*n]);
    u3 = x([8*n+1:9*n]);

    for i=[1:n-1]
        obj = obj + 1/2 * h * (u1(i)^2 + u1(i+1)^2 + u2(i)^2 + u2(i+1)^2 + u3(i)^2 + u3(i+1)^2);
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
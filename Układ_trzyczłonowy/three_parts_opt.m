clear;
clc;

# Parametry modelu:
# Masy, długości i momenty bezwładności członów
global l1=1;
global m1=1;
global J1=1/12 * m1 * l1^2;
global l2=1;
global m2=1;
global J2=1/12 * m2 * l2^2;
global l3=1;
global m3=1;
global J3=1/12 * m3 * l3^2;

# Przyśpieszenie grawitacyjne
global grav=0;

# Położenie początkowe i końcowe końcówki ostatniego członu
global xp=2;
global yp=0.5;
global xk=2.5;
global yk=1;

# Czas symulacji
T=1;

# Liczba kroków czasowych
global n=16;

# Maksymalna liczba iteracji algorytmu SQP
max_iter=300;


tic();

function rot=R(fi)
    rot=[cos(fi), -sin(fi); sin(fi), cos(fi)];
endfunction

global h=T/(n-1);
t=[0:h:T];

#Stałe:
global I=[1, 0; 0, 1];
global Q=[0, -1; 1, 0]; 

global M=diag([m1;m1;J1;m2;m2;J2;m3;m3;J3]);
#global s_A1=[-l1/2;0];
#global s_B1=[l1/2;0];
#global s_B2=[-l2/2;0];
#global s_C2=[l2/2;0];
#global s_C3=[-l3/2;0];
#global s_D3=[l3/2;0];


l=[l1,l2,l3];

global S1=zeros(2,3);
global S2=zeros(2,3);
for part=[1:1:3]
      S1(:,part)=[-l(part)/2;0];
      S2(:,part)=[l(part)/2;0];
endfor

#Dla bez g:
theta10=0.4151+t'/T*(0.4037-0.4151);
theta20=0.6687+t'/T*(0.5246-(0.6687));
theta30=-1.9897+t'/T*(-1.1228-(-1.9897));
x0=[theta10;theta20;theta30;zeros(3*n,1);ones(3*n,1)];
#Z g:
#theta10=1-t'/T;
#theta20=-1-t'/T;
#theta30=-1+t'/T;

#x0=[zeros(6*n,1);ones(3*n,1)];
#x0=[zeros(9*n,1);];
#x0=[theta10;theta20;theta30;zeros(3*n,1);20*ones(n,1);5*ones(n,1);-5*ones(n,1)];
#load zg4.mat
#x0([2*n+1:3*n])=-1*ones(n,1);
#x0([8*n+1:9*n])=zeros(n,1);

load tvp2.mat
k = rows(x)/9;

tp=[0:T/(k-1):T]';

p10 = interp1(tp, x([1:k]), t');
p20 = interp1(tp, x([k+1:2*k]), t');
p30 = interp1(tp, x([2*k+1:3*k]), t');
v10 = interp1(tp, x([3*k+1:4*k]), t');
v20 = interp1(tp, x([4*k+1:5*k]), t');
v30 = interp1(tp, x([5*k+1:6*k]), t');
u10 = interp1(tp, x([6*k+1:7*k]), t');
u20 = interp1(tp, x([7*k+1:8*k]), t');
u30 = interp1(tp, x([8*k+1:9*k]), t');
clear x
#x0 = [p10; p20; p30; v10; v20; v30; u10; u20; u30];
#x0=[zeros(9*n,1);];

function r = f_ogr(x)
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

    #r =[thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n); thetadot3(1); thetadot3(n); l1*cos(theta1(n))+l2*cos(theta2(n)+theta1(n))+l3*cos(theta3(n)+theta2(n)+theta1(n))-xk; l1*sin(theta1(n))+l2*sin(theta2(n)+theta1(n))+l3*sin(theta3(n)+theta2(n)+theta1(n))-yk; l1*cos(theta1(1))+l2*cos(theta2(1)+theta1(1))+l3*cos(theta3(1)+theta2(1)+theta1(1))-xp; l1*sin(theta1(1))+l2*sin(theta2(1)+theta1(1))+l3*sin(theta3(1)+theta2(1)+theta1(1))-yp]; # Warunki brzegowe

    r=[theta1(1)-0.4151; theta1(n)-0.4037; theta2(1)-0.6687; theta2(n)-0.5246; theta3(1)+1.9897; theta3(n)+1.1228;
    thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n); thetadot3(1); thetadot3(n);];

    for i=[1:n-1]
        vec = q(:,i+1)-q(:,i)-1/2*h*(three_parts_dynamics_function(q(:,i+1),u(:,i+1))+three_parts_dynamics_function(q(:,i),u(:,i)));
        r=[r;vec];
    endfor

endfunction

function r = bun(x)
    global n;
    theta1 = x([1:n]);
    theta2 = x([n+1:2*n]);
    theta3 = x([2*n+1:3*n]);

    r=[theta1+pi;-theta1+pi;theta2+pi;-theta2+pi;theta3+pi;-theta3+pi];

endfunction

function obj = f_celu(x)
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
#lb=[-100*ones(3*n,1);-100*ones(6*n,1)];
#ub=[100*ones(3*n,1);100*ones(6*n,1)];
#[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, @bun,[],[],max_iter);
[x, obj, info, iter, nf, lambda] = sqp (x0, @f_celu, @f_ogr, [],[],[],max_iter);
time_of_execution=toc()

obj
info
iter

#save n5.mat x

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
#plot(l1*cos(theta1)+l2*cos(theta2+theta1)+l3*cos(theta3+theta2+theta1), l1*sin(theta1)+l2*sin(theta2+theta1)+l3*sin(theta3+theta2+theta1), "LineWidth",3)
#axis("equal")
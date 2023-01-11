clear;
clc;

function rot=R(fi)
    rot=[cos(fi), -sin(fi); sin(fi), cos(fi)];
endfunction

T=4;
h=0.001;
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

global M=diag([m1;m1;J1;m2;m2;J2]);
global s_A1=[-l1/2;0];
global s_B1=[l1/2;0];
global s_B2=[-l2/2;0];
global s_C2=[l2/2;0];

l=[l1,l2];

global S1=zeros(2,2);
global S2=zeros(2,2);
for part=[1:1:2]
      S1(:,part)=[-l(part)/2;0];
      S2(:,part)=[l(part)/2;0];
endfor


global grav=-9.81;

theta1=0;
theta2=0;
thetadot1=0;
thetadot2=0;

x0=[theta1,theta2,thetadot1,thetadot2];

function r=d(part, joint, S1_g, S2_g)
   global Q
   r=[0;0];
   for i=[joint:1:part]
      r=r-S1_g(:,i)+S2_g(:,i);
   endfor
   r=r-S2_g(:,part);
   r=Q*r;
endfunction

function r=d_dot(part, joint, S1_g, S2_g, fidot)
   global Q
   r=[0;0];
   for i=[joint:1:part]
      r=r+S1_g(:,i)*fidot(i)-S2_g(:,i)*fidot(i);
   endfor
   r=r+S2_g(:,part)*fidot(part);
endfunction

function r=fcn(x,t)
    global I
    global Q 
    global M
    global grav
    global m1
    global m2
    global S1
    global S2

    #t

    g=[0;grav*m1;0;0;grav*m2;0];

    #g=[0;grav*m1;mom1(t)-mom2(t);0;grav*m2;mom2(t)];

    #[0;grav*m1;mom1-mom2;0;grav*m2;mom2]
 
    thetadot=[x(3);x(4)];

    fi=[x(1);x(2)+x(1)];
    fidot=[x(3);x(3)+x(4)];

    S1_g=zeros(2,2);
    S2_g=zeros(2,2);
    for part=[1:1:2]
      S1_g(:,part)=R(fi(part))*S1(:,part);
      S2_g(:,part)=R(fi(part))*S2(:,part);
   endfor

   D=[I, Q * S1_g(:,1), zeros(2,3);
      I, Q * S2_g(:,1), -I, -Q * S1_g(:,2)];


   B=[d(1,1,S1_g,S2_g), zeros(2,1);
        1, 0;
       d(2,1,S1_g,S2_g), d(2,2,S1_g,S2_g);
          1, 1];

    #D*B   

   Bdot=[d_dot(1,1,S1_g,S2_g,fidot), zeros(2,1);
        0, 0;
       d_dot(2,1,S1_g,S2_g,fidot), d_dot(2,2,S1_g,S2_g,fidot);
          0, 0];


    M_new=B' * M * B;
    g_new=B' * (g - M * Bdot * thetadot);

    r1=x(3);
    r2=x(4);
    r3=M_new\g_new;
    r=[r1,r2,r3'];

endfunction

[x, istate, msg] = lsode (@fcn, x0, t);

g=[0;grav*m1;1;0;grav*m2;1];

Ek=zeros(1,n);
Ep=zeros(1,n);
Ec=zeros(1,n);
x1=zeros(1,n);
y1=zeros(1,n);
xb=zeros(1,n);
yb=zeros(1,n);
xc=zeros(1,n);
yc=zeros(1,n);

for i=[1:1:n]
    fi1=x(i,1);
    fi2=x(i,2)+x(i,1);
    fidot1=x(i,3);
    fidot2=x(i,4)+x(i,3);  
    thetadot=[x(i,3);x(i,4)];

    s_A1_g=R(fi1)*s_A1;
    s_B1_g=R(fi1)*s_B1;
    s_B2_g=R(fi2)*s_B2;
    s_C2_g=R(fi2)*s_C2;

    D=[I, Q * s_A1_g, zeros(2,3);
       I, Q * s_B1_g, -I, -Q * s_B2_g];
    

    d_11=-Q * s_A1_g;
    d_21=Q*(-s_A1_g+s_B1_g-s_B2_g);
    d_22=-Q*s_B2_g;

    B=[d_11, zeros(2,1);
          1, 0;
       d_21, d_22;
          1, 1];
    
    v=B*thetadot;
    q=[-s_A1_g;fi1;-s_A1_g+s_B1_g-s_B2_g;fi2];
    x1(i)=-s_A1_g(1);
    y1(i)=-s_A1_g(2);
    xb(i)=-s_A1_g(1)+s_B1_g(1);
    yb(i)=-s_A1_g(2)+s_B1_g(2);
    xc(i)=-s_A1_g(1)+s_B1_g(1)-s_B2_g(1)+s_C2_g(1);
    yc(i)=-s_A1_g(2)+s_B1_g(2)-s_B2_g(2)+s_C2_g(2);

    Ek(i)=1/2*v'*M*v;
    Ep(i)=-q'*g;
    Ec(i)=Ek(i)+Ep(i);
endfor
plot(t,x1,t,y1)
#plot(xc,yc,xks,yks)
#plot([zeros(1,101);xb],[zeros(1,101);yb],[xb;xc],[yb;yc])

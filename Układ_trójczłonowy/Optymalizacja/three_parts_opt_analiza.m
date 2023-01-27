clear;
clc;

### Parametry modelu: ###
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

### Koniec parametrów ###

# Dyskretyzacja czasu
global h=T/(n-1);
t=[0:h:T];

# Macierz masowa
global M=diag([m1;m1;J1;m2;m2;J2;m3;m3;J3]);

# Wektor długości
l=[l1,l2,l3];

# Wektory od środka masy do złącz
global P=zeros(2,3);
global S=zeros(2,3);
for part=[1:1:3]
      P(:,part)=[-l(part)/2;0];
      S(:,part)=[l(part)/2;0];
endfor


function r = f_ogr(x)
    global n;
    global h;
    global l1;
    global l2;
    global l3;
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

    for i=[1:n-1]
        vec = q(:,i+1)-q(:,i)-1/2*h*(three_parts_dynamics_function(q(:,i+1),u(:,i+1))+three_parts_dynamics_function(q(:,i),u(:,i)));
        r=[r;vec];
    endfor

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

T=1;
exe=[];
objs=[];
iters=[];
for j=[3:1:40]
    clear global n
    clear global h
    global n=j
    global h=T/(n-1);
    t=[0:h:T];
    j
    n
    theta10=0.4151+t'/T*(-0.4037-0.4151);
    theta20=0.6687+t'/T*(0.5246-(0.6687));
    theta30=-1.9897+t'/T*(-1.1228-(-1.9897));

    x0=[theta10;theta20;theta30;zeros(3*n,1);ones(3*n,1)];

    tic();
    [x, obj, info, iter, nf, lambda] = sqp (x0, @f_celu, @f_ogr, [],[],[],max_iter);
    time_of_execution=toc();
    exe=[exe,time_of_execution];
    objs=[objs,obj];
    iters=[iters, iter];
endfor
j=[3:1:40];
plot(j,exe,j,objs,j,iters)
save analiza.mat j exe objs iters
save n40.mat x

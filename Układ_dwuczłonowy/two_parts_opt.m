clear;
clc;

### Parametry modelu: ###
# Masy, długości i momenty bezwładności członów
global l1=1;
global m1=1;
global J1=1;
#global J1=1/12 * m1 * l1^2;
global l2=1;
global m2=1;
global J2=1;
#global J2=1/12 * m2 * l2^2;

# Przyśpieszenie grawitacyjne
global grav=0;

# Położenie początkowe i końcowe końcówki ostatniego członu
#global xp=2;
#global yp=0.5;
#global xk=2.5;
#global yk=1;

# Położenie początkowe i końcowe członów (współrzędne złączowe)
global theta1p=0;
global theta2p=0;
global theta1k=pi/2;
global theta2k=pi/2;

# Czas symulacji
T=1;

# Liczba kroków czasowych
global n=21;

# Maksymalna liczba iteracji algorytmu SQP
max_iter=300;

### Koniec parametrów ###

# Dyskretyzacja czasu
global h=T/(n-1);
t=[0:h:T];

# Macierz masowa
global M=diag([m1;m1;J1;m2;m2;J2]);

# Wektor długości
l=[l1,l2];

# Wektory od środka masy do złącz
global P=zeros(2,2);
global S=zeros(2,2);
for part=[1:1:2]
      P(:,part)=[-l(part)/2;0];
      S(:,part)=[l(part)/2;0];
endfor

# Przybliżenie początkowe
x0=ones(6*n,1);

function r = f_ogr(x)
    global n;
    global h;
    global l1;
    global l2;
    #global xk;
    #global yk;
    #global xp;
    #global yp;
    global theta1p;
    global theta2p;
    global theta1k;
    global theta2k;

    theta1 = x([1:n]);
    theta2 = x([n+1:2*n]);
    thetadot1 = x([2*n+1:3*n]);
    thetadot2 = x([3*n+1:4*n]);
    q=[theta1';theta2';thetadot1';thetadot2'];
    u=[x([4*n+1:5*n])';x([5*n+1:6*n])'];

    #r =[thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n); l1*cos(theta1(n))+l2*cos(theta2(n)+theta1(n))-xk; l1*sin(theta1(n))+l2*sin(theta2(n)+theta1(n))-yk; l1*cos(theta1(1))+l2*cos(theta2(1)+theta1(1))-xp; l1*sin(theta1(1))+l2*sin(theta2(1)+theta1(1))-yp]; # Warunki brzegowe
    r =[theta1(1)-theta1p; theta1(n)-theta1k; theta2(1)-theta2p; theta2(n)-theta2k; thetadot1(1); thetadot1(n); thetadot2(1); thetadot2(n)];
    for i=[1:n-1]
        vec = q(:,i+1)-q(:,i)-1/2*h*(two_parts_dynamics_function(q(:,i+1),u(:,i+1))+two_parts_dynamics_function(q(:,i),u(:,i)));
        r=[r;vec];
    endfor

endfunction

function obj = f_celu(x)
    global n;
    global h;
    obj = 0;
    u1 = x([4*n+1:5*n]);
    u2 = x([5*n+1:6*n]);

    for i=[1:n-1]
        obj = obj + 1/2 * h * (u1(i)^2 + u1(i+1)^2 + u2(i)^2 + u2(i+1)^2);
    endfor

endfunction

tic();
[x, obj, info, iter, nf, lambda] = sqp (x0, @f_celu, @f_ogr, [],[],[],max_iter);
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

#plot(t,theta1,t,theta2)

save two_parts_results.mat x

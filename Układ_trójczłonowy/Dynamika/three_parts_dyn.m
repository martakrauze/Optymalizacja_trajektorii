clear;
clc;

# Czas symulacji
T=10;
h=0.01;
t=[0:h:T];
n=length(t)

# Parametry układu
global l1=1;
global m1=1;
global J1=1/12 * m1 * l1^2;
global l2=1;
global m2=1;
global J2=1/12 * m2 * l2^2;
global l3=1;
global m3=1;
global J3=1/12 * m3 * l3^2;

# Grawitacja
global grav=-9.81;

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


# Stan początkowy
theta1=0;
theta2=0;
theta3=0;
thetadot1=0;
thetadot2=0;
thetadot3=0;

x0=[theta1,theta2, theta3, thetadot1,thetadot2, thetadot3];

# Rozwiązywanie układu równań różniczkowych
[x, istate, msg] = lsode (@three_parts_dynamics_function, x0, t);


fi1=x(:,1);
fi2=x(:,2)+x(:,1);
fi3=x(:,3)+x(:,2)+x(:,1);
fidot1=x(:,4);
fidot2=x(:,4)+x(:,5);
fidot3=x(:,4)+x(:,5)+x(:,6);
x3=l1*cos(fi1)+l2*cos(fi2)+l3/2*cos(fi3);
y3=l1*sin(fi1)+l2*sin(fi2)+l3/2*sin(fi3);

# Wczytanie przygotowanych danych otrzymanych za pomocą programu FreeDyn
load t.txt
load x3.txt
load y3.txt

plot(t(1:1:401),x3(1:1:401),"LineWidth",3,t(1:1:401),y3(1:1:401),"LineWidth",3,t_free(1:1:401),x_free(1:1:401),"LineWidth",2,t_free(1:1:401),y_free(1:1:401),"LineWidth",2)
legend("3.part x - Octave", "3.part y - Octave","3.part x - FreeDyn","3.part y - FreeDyn")
xlabel("czas[s]")
ylabel("położenie[m]")
set(gca, "fontsize", 20)


#plot(t(1:1:401),abs(x3(1:1:401)-x_free(1:1:401)),"LineWidth",3,t(1:1:401),abs(y3(1:1:401)-y_free(1:1:401)),"LineWidth",3)
#set(gca, "fontsize", 20)
#legend("moduł różnicy x", "moduł różnicy y","location","northwest")
#xlabel("czas[s]")

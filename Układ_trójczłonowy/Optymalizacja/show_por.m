clear;
clc;

load n5.mat
T=1;
global n=length(x)/9;
global h=T/(n-1);
t=[0:h:T];
load n10.mat
global n2=length(x2)/9;
global h2=T/(n2-1);
t2=[0:h2:T];

load n40.mat
global n3=length(x3)/9;
global h3=T/(n3-1);
t3=[0:h3:T];

global l1=1;
global l2=1;
global l3=1;

N=3

l=[l1,l2,l3];
#l=[l1,l2];

p=[]; v=[]; u=[];

for i=[1:N]
    b = x([(i-1)*n+1:i*n]);
    p=[p,b];
    b = x([n*N+(i-1)*n+1:n*N+i*n]);
    v=[v,b];
    b = x([2*n*N+(i-1)*n+1:2*n*N+i*n]);
    u=[u,b];
endfor

p2=[]; v2=[]; u2=[];

for i=[1:N]
    b = x2([(i-1)*n2+1:i*n2]);
    p2=[p2,b];
    b = x2([n2*N+(i-1)*n2+1:n2*N+i*n2]);
    v2=[v2,b];
    b = x2([2*n2*N+(i-1)*n2+1:2*n2*N+i*n2]);
    u2=[u2,b];
endfor

p3=[]; v3=[]; u3=[];

for i=[1:N]
    b = x3([(i-1)*n3+1:i*n3]);
    p3=[p3,b];
    b = x3([n3*N+(i-1)*n3+1:n3*N+i*n3]);
    v3=[v3,b];
    b = x3([2*n3*N+(i-1)*n3+1:2*n3*N+i*n3]);
    u3=[u3,b];
endfor
#p=p+ones(n,N)*6*pi

fi=[zeros(n,N)];
for part=[1:N]
    for i=[1:n]
        for j=[1:part]
            fi(i,part)=fi(i,part)+p(i,j);
        endfor
    endfor
endfor

trajx = [zeros(1,n)];
trajy = [zeros(1,n)];
for i=[1:N]
    trajx = [trajx; trajx(i,:) + (l(i) * cos(fi(:,i)))'];
    trajy = [trajy; trajy(i,:) + (l(i) * sin(fi(:,i)))'];
endfor

xkon=l*cos(fi)';
ykon=l*sin(fi)';

fi2=[zeros(n2,N)];
for part=[1:N]
    for i=[1:n2]
        for j=[1:part]
            fi2(i,part)=fi2(i,part)+p2(i,j);
        endfor
    endfor
endfor

trajx2 = [zeros(1,n2)];
trajy2 = [zeros(1,n2)];
for i=[1:N]
    trajx2 = [trajx2; trajx2(i,:) + (l(i) * cos(fi2(:,i)))'];
    trajy2 = [trajy2; trajy2(i,:) + (l(i) * sin(fi2(:,i)))'];
endfor

xkon2=l*cos(fi2)';
ykon2=l*sin(fi2)';

fi3=[zeros(n3,N)];
for part=[1:N]
    for i=[1:n3]
        for j=[1:part]
            fi3(i,part)=fi3(i,part)+p3(i,j);
        endfor
    endfor
endfor

trajx3 = [zeros(1,n3)];
trajy3 = [zeros(1,n3)];
for i=[1:N]
    trajx3 = [trajx3; trajx3(i,:) + (l(i) * cos(fi3(:,i)))'];
    trajy3 = [trajy3; trajy3(i,:) + (l(i) * sin(fi3(:,i)))'];
endfor

xkon3=l*cos(fi3)';
ykon3=l*sin(fi3)';

function plotangle(t,p,t2,p2)
plot(t,p,"LineWidth", 2,"black", t2, p2,"LineWidth", 2,"blue");
title("kąt")
xlabel("t[s]")
ylabel('\theta [rad]', 'interpreter', 'tex')
legend("człon1","człon2", "człon3")
set(gca, "fontsize", 17)
grid
endfunction

function plotvel(t,v,t2,v2)
#subplot(3,2,3)
plot(t,v,"LineWidth", 2, "black", t2, v2,"LineWidth", 2,"blue");
title("prędkość kątowa")
xlabel("t[s]")
ylabel('\omega [rad/s]', 'interpreter', 'tex')
legend("człon1","człon2", "człon3")
set(gca, "fontsize", 17)
grid
endfunction

function plotmom(t,u,t2,u2)
#subplot(3,2,5)
plot(t,u,"LineWidth", 2, "black", t2, u2,"LineWidth", 2,"blue");
title("moment zadany")
xlabel("t[s]")
ylabel("M[Nm]")
legend("człon1","człon2", "człon3")
set(gca, "fontsize", 17)
grid
endfunction

function plottraj(trajx,trajy,xkon,ykon, trajx2,trajy2,xkon2,ykon2, trajx3,trajy3,xkon3,ykon3)
global n
#subplot(3,2,2:2:6)
#plot(trajx(:,1),trajy(:,1),"LineWidth", 4,trajx(:,n),trajy(:,n),"black","LineWidth", 4)
#plot(trajx,trajy,"black","LineWidth", 4)
#
plot(xkon3,ykon3,"LineWidth", 2,"black")
hold;
plot(xkon2,ykon2,"LineWidth", 2,"red")
plot(xkon,ykon,"LineWidth", 2,"blue")
plot(xkon3,ykon3,"o", "markeredgecolor", "black", "markerfacecolor", "black")
plot(xkon2,ykon2,"o", "markeredgecolor", "red", "markerfacecolor", "red")
plot(xkon,ykon,"o", "markeredgecolor", "blue", "markerfacecolor", "blue")
title("trajektoria")
legend("n=40","n=10","n=5")
xlabel("x[m]")
ylabel("y[m]")
axis("equal")
set(gca, "fontsize", 17)
grid
endfunction

#subplot(3,1,1)
#plotangle(t,p,t2,p2)
#subplot(3,1,2)
#plotvel(t,v,t2,v2)
#subplot(3,1,3)
#plotmom(t,u,t2,u2)
plottraj(trajx,trajy,xkon,ykon, trajx2,trajy2,xkon2,ykon2, trajx3,trajy3,xkon3,ykon3)
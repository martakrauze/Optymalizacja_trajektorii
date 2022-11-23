clear;
clc;

load wyniki.mat

n = rows(x)/(4*N);

p=[]; v=[]; a=[]; u=[];

for i=[1:N]
    b = x([(i-1)*n+1:i*n]);
    p=[p,b];
    b = x([n*N+(i-1)*n+1:n*N+i*n]);
    v=[v,b];
    b = x([2*n*N+(i-1)*n+1:2*n*N+i*n]);
    a=[a,b];
    b = x([3*n*N+(i-1)*n+1:3*n*N+i*n]);
    u=[u,b];
endfor

trajx = [zeros(1,n)];
trajy = [zeros(1,n)];
for i=[1:N]
    trajx = [trajx; trajx(i,:) + (l(i) * cos(p(:,i)))'];
    trajy = [trajy; trajy(i,:) + (l(i) * sin(p(:,i)))'];
endfor

xkon=l*cos(p)';
ykon=l*sin(p)';

subplot(4,2,1)
plot(t,p);
title("kąt")
xlabel("t[s]")
ylabel("fi[rad]")
legend("człon1","człon2", "człon3")

subplot(4,2,3)
plot(t,v);
title("prędkość kątowa")
xlabel("t[s]")
ylabel("omega[rad/s]")
legend("człon1","człon2", "człon3")

subplot(4,2,5)
plot(t,a);
title("przyspieszenie kątowe")
xlabel("t[s]")
ylabel("epsilon[rad/s^2]")
legend("człon1","człon2", "człon3")

subplot(4,2,7)
plot(t,u);
title("moment zadany")
xlabel("t[s]")
ylabel("M[Nm]")
legend("człon1","człon2", "człon3")

subplot(4,2,2:2:8)
plot(trajx,trajy)
hold;
plot(xkon,ykon,"*")
title("trajektoria")
xlabel("x[m]")
ylabel("y[m]")
axis("equal")
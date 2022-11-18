clear;
clc;

load wyniki.mat

n = rows(x)/8;

# Wyświetlanie wyników 

p1 = x([1:n]);
p2 = x([n+1:2*n]);
v1 = x([2*n+1:3*n]);
v2 = x([3*n+1:4*n]);
a1 = x([4*n+1:5*n]);
a2 = x([5*n+1:6*n]);
u1 = x([6*n+1:7*n]);
u2 = x([7*n+1:8*n]);

trajx = [zeros(1,n); (l1 * cos(p1))'; (l1 * cos(p1))' + (l2 * cos(p2))'];
trajy = [zeros(1,n); (l1 * sin(p1))'; (l1 * sin(p1))' + (l2 * sin(p2))'];


xkon=l1*cos(p1)+l2*cos(p2);
ykon=l1*sin(p1)+l2*sin(p2);

#subplot(4,2,2:2:8)
plot(trajx,trajy)
title("trajektoria")
xlabel("x[m]")
ylabel("y[m]")
axis("equal")
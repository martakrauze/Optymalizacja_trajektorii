#Block-Move Example na podstawie publikacji M. Kelley https://epubs.siam.org/doi/10.1137/16M1062569

clear;
clc;

# krok całkowania
global n=40;
global h=1/(n-1);

# p - przemieszczenie, v - prędkość, a - przyspieszenie

# przybliżenie początkowe - przmieszczenie liniowo, prędkość = 1, przyspieszenie = 0
p0 = zeros(n,1);
for i=[1:n]
    p0(i)=h*(i-1);
endfor
v0 = ones(n,1);
a0 = zeros(n,1);
x0 = [p0; v0; a0];

function vec = trapezoidal_transcription(x, xdot, n, h)

    vec = zeros(n-1,1);
    for i=[1:n-1]
        vec(i) = x(i+1)-x(i)-1/2*h*(xdot(i+1)+xdot(i));
    endfor

endfunction

# Zdyskretyzowanie równania ruchu(równowagi)
function r = g (x)

    global n;
    global h;
    p = x([1:n]);
    v = x([n+1:2*n]);
    a = x([2*n+1:3*n]);

    r0 =[p(1); p(n)-1; v(1); v(n)]; # Warunki brzegowe

    r1 = trapezoidal_transcription(p, v, n, h);
    r2 = trapezoidal_transcription(v, a, n, h);

    r = [r0; r1; r2];

endfunction

# Funkcja celu - całka po czasie z kwadratu przypieszenia
# przyspieszenie = siła / masa - zatem minimalizuje siłę

function obj = phi (x)

    global n;
    global h;
    obj = 0;
    a = x([2*n+1:3*n]);
    
    for i=[1:n-1]
        obj = obj + 1/2 * h * (a(i)^2 + a(i+1)^2);
    endfor

endfunction

# Algorytm optymalizacyjny
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, []);

obj
info
iter

# Wyświetlanie i porównanie wyników z teoretycznymi
t=[0:h:1];

p = x([1:n]);
v = x([n+1:2*n]);
a = x([2*n+1:3*n]);

pstar=(3*t.^2-2*t.^3)';
vstar=(6*t-6*t.^2)';
astar=(6-12*t)';

subplot(3,1,1)
plot(t,p,t,pstar);
title("przemieszczenie")
xlabel("t[s]")
ylabel("x[m]")
legend("numeryczne","analityczne")

subplot(3,1,2)
plot(t,v,t,vstar);
title("prędkość")
xlabel("t[s]")
ylabel("v[m/s]")
legend("numeryczne","analityczne")

subplot(3,1,3)
plot(t,a,t,astar);
title("przyspieszenie")
xlabel("t[s]")
ylabel("a[m/s^2]")
legend("numeryczne","analityczne")
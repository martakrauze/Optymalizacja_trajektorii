#Cart-Pole Swing-Up Example na podstawie publikacji M. Kelley https://epubs.siam.org/doi/10.1137/16M1062569

clear;
clc;

# krok całkowania
global n = 25;
global T = 2; #czas wykonania manewru
global h = T/(n-1);

#zmienne globalne
global dist = 1; #miejsce wykonania manewru
global m1 = 1; #masa wózka
global m2 = 0.3; #masa wachadła
global grav = 9.81; #przyspieszenie ziemskie
global len = 0.5; #długość wachadła
#umax=20
#dmax=2

# p - przemieszczenie, v - prędkość, a - przyspieszenie, 1- wózek, 2- wachadło

# przybliżenie początkowe - przmieszczenie liniowo, prędkość = 1, przyspieszenie = 0

q10 = zeros(n,1);
q20 = zeros(n,1);
v10 = zeros(n,1);
v20 = zeros(n,1);
a10 = zeros(n,1);
a20 = zeros(n,1);
u0 = zeros(n,1);

for i=[1:n]
    q10(i)=h*(i-1)/T * dist;
endfor
for i=[1:n]
    q20(i)=h*(i-1)/T * pi;
endfor
x0 = [q10;q20;v10;v20;a10;a20;u0];

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
    global dist;
    global m1;
    global m2;
    global grav;
    global len;
   
    q1 = x([1:n]);
    q2 = x([n+1:2*n]);
    v1 = x([2*n+1:3*n]);
    v2 = x([3*n+1:4*n]);
    a1 = x([4*n+1:5*n]);
    a2 = x([5*n+1:6*n]);
    u = x([6*n+1:7*n]);

    r0= [q1(1); q2(1); q1(n)-dist; q2(n)-pi; v1(1); v2(1); v1(n); v2(n)]; # Warunki brzegowe

    r1 = trapezoidal_transcription(q1, v1, n, h);
    r2 = trapezoidal_transcription(q2, v2, n, h);
    r3 = trapezoidal_transcription(v1, a1, n, h);
    r4 = trapezoidal_transcription(v2, a2, n, h);
    r5 = (len * m2 * sin(q2) .* v2.^2 + u + m2 * grav * cos(q2) .* sin(q2)) ./ (m1 + m2 * (1 - cos(q2).^2)) - a1;
    r6 = - (len * m2 * cos(q2) .* sin(q2) .* v2.^2 + u .* cos(q2) + (m1 + m2) * grav * sin(q2) ) ./ (len * m1 + len * m2 * (1 - cos(q2).^2)) - a2;

    r=[r0; r1; r2; r3; r4; r5; r6];

endfunction

# Funkcja celu - całka po czasie z kwadratu siły
function obj = phi (x)

    global n;
    global h;
    obj = 0;
    u = x([6*n+1:7*n]);
    for i=[1:n-1]
        obj = obj + 1/2 * h * (u(i)^2 + u(i+1)^2);
    endfor

endfunction

# Algorytm optymalizacyjny
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, []);

obj
info
iter

# Wyświetlanie wyników 

t=[0:h:T];

q1 = x([1:n]);
q2 = x([n+1:2*n]);
v1 = x([2*n+1:3*n]);
v2 = x([3*n+1:4*n]);
a1 = x([4*n+1:5*n]);
a2 = x([5*n+1:6*n]);
u = x([6*n+1:7*n]);

trajx = [ q1'; (q1 + sin(q2))'];
trajy = [zeros(1,n); (-cos(q2))'];


subplot(4,2,1)
plot(t,q1);
title("przemieszczenie")
xlabel("t[s]")
ylabel("x[m]")

subplot(4,2,3)
plot(t,q2);
title("kąt")
xlabel("t[s]")
ylabel("fi[rad]")

subplot(4,2,5)
plot(t,u);
title("siła")
xlabel("t[s]")
ylabel("F[N]")

subplot(4,2,2:2:6)
plot(trajx,trajy)
title("trajektoria")
xlabel("x[m]")
ylabel("y[m]")
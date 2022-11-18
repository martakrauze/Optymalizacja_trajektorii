clear;
clc;

# krok całkowania
global n=10;
global T=1;
global h=T/(n-1);
global m1=1; #masa
global l1=1; #długość
global m2=1; #masa
global l2=1.2; #długość
# global grav=9.81;

# położenie końcowe końcówki układu

global xk=0.5;
global yk=1.5;

# p - kąt, v - prędkość kątowa, a - przyspieszenie kątowe, u - moment napędzający

# załadowanie poprzedniego wyniku lub przybliżenia początkowego

load wyniki.mat

k = rows(x)/8;

tp=[0:T/(k-1):T]';

t=[0:h:T]';

p10 = interp1(tp, x([1:k]), t);
p20 = interp1(tp, x([k+1:2*k]), t);
v10 = interp1(tp, x([2*k+1:3*k]), t);
v20 = interp1(tp, x([3*k+1:4*k]), t);
a10 = interp1(tp, x([4*k+1:5*k]), t);
a20 = interp1(tp, x([5*k+1:6*k]), t);
u10 = interp1(tp, x([6*k+1:7*k]), t);
u20 = interp1(tp, x([7*k+1:8*k]), t);

clear x;

x0 = [p10; p20; v10; v20; a10; a20; u10; u20];

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
    global m1;
    global m2;
    global l1;
    global l2;
    global grav;
    global xk;
    global yk;

    p1 = x([1:n]);
    p2 = x([n+1:2*n]);
    v1 = x([2*n+1:3*n]);
    v2 = x([3*n+1:4*n]);
    a1 = x([4*n+1:5*n]);
    a2 = x([5*n+1:6*n]);
    u1 = x([6*n+1:7*n]);
    u2 = x([7*n+1:8*n]);

    r0 =[p1(1); p2(1); v1(1); v1(n); v2(1); v2(n); l1*cos(p1(n))+l2*cos(p2(n))-xk; l1*sin(p1(n))+l2*sin(p2(n))-yk]; # Warunki brzegowe

    r1 = trapezoidal_transcription(p1, v1, n, h);
    r2 = trapezoidal_transcription(p2, v2, n, h);
    r3 = trapezoidal_transcription(v1, a1, n, h);
    r4 = trapezoidal_transcription(v2, a2, n, h);
    
    #brak grawitacji
    r5 = 1/3 * m2 * l2^2 * a2 - u2; #+ l2/2 * cos(p2) * m2 * grav - u2;
    r6 = 1/3 * m1 * l1^2 * a1 + (m2 * l2 * l1)/2 * (sin(p1) .* cos(p2) .* v2.^2 + sin(p1) .* sin(p2) .* a2
    - cos(p1) .* sin(p2) .* v2.^2 + cos(p1) .* cos(p2) .* a2) - u1;

    r = [r0; r1; r2; r3; r4; r5; r6];

endfunction

# Funkcja celu - całka po czasie z sumy kwadratów momentów
function obj = phi (x)
    global n;
    global h;
    obj = 0;
    u1 = x([6*n+1:7*n]);
    u2 = x([7*n+1:8*n]);

    for i=[1:n-1]
        obj = obj + 1/2 * h * (u1(i)^2 + u1(i+1)^2 + u2(i)^2 + u2(i+1)^2);
    endfor

endfunction

# Algorytm optymalizacyjny
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, [],[],[],300);

obj
info
iter

save wyniki.mat x t l1 l2 #zapis wyników

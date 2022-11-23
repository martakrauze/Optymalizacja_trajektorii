clear;
clc;

#liczba członów
global N=3;

# krok całkowania
global n=5;
global T=1; #czas ruchu
global h=T/(n-1);
global m=[1,1,1]; #masa
global l=[1,1.2,0.5]; #długość
# global grav=9.81;

# położenie końcowe końcówki układu

global xk=1;
global yk=1;

# p - kąt, v - prędkość kątowa, a - przyspieszenie kątowe, u - moment napędzający

x0 = zeros(4*N*n,1);

function vec = trapezoidal_transcription(x, xdot, n, h)
    vec = zeros(n-1,1);
    for i=[1:n-1]
        vec(i) = x(i+1)-x(i)-1/2*h*(xdot(i+1)+xdot(i));
    endfor
endfunction


# Zdyskretyzowanie równania ruchu(równowagi)
function r = g (x)
    global n;
    global N;
    global h;
    global m;
    global l;
    global xk;
    global yk;

    p=[]; v=[]; a=[]; u=[];

    w=[]; r0=[];

    for i=[1:N]
        b = x([(i-1)*n+1:i*n]);
        p=[p,b];
        b = x([n*N+(i-1)*n+1:n*N+i*n]);
        v=[v,b];
        b = x([2*n*N+(i-1)*n+1:2*n*N+i*n]);
        a=[a,b];
        b = x([3*n*N+(i-1)*n+1:3*n*N+i*n]);
        u=[u,b];

        b = trapezoidal_transcription(p(:,i), v(:,i), n, h);
        w = [w;b];
        b = trapezoidal_transcription(v(:,i), a(:,i), n, h);
        w = [w;b];
        r0 = [r0; p(1,i); v(1,i); v(n,i);];
    endfor

    r0 =[r0; l*cos(p(n,:))'-xk; l*sin(p(n,:))'-yk]; # Warunki brzegowe
    
    r1=[];

    for j=[1:n]
    for i=[N:-1:1]
        if i<N
        S(:,i)=F(:,i+1)+S(:,i+1);
        else
        S(:,i)=zeros(2,1);
        endif
        F(1,i)=-l(i)/2 * (cos(p(j,i)) .* v(j,i).^2 + sin(p(j,i)) .* a(j,i));
        F(2,i)=l(i)/2 * (-sin(p(j,i)) .* v(j,i).^2 + cos(p(j,i)) .* a(j,i));
        vec=[-l(i)*sin(p(j,i)), l(i)*cos(p(j,i))];
        b=1/3 * m(i) * l(i)^2 * a(j,i) - u(j,i) + vec * S(:,i);
        r1=[r1;b];
    endfor
    endfor

    r = [r0; w; r1];

endfunction

# Funkcja celu - całka po czasie z sumy kwadratów momentów
function obj = phi (x)
    global n;
    global N;
    global h;
    obj = 0;

    for i=[1:N]
        u = x([3*n*N+(i-1)*n+1:3*n*N+i*n]);
        for j=[1:n-1]
            obj = obj + 1/2 * h * (u(j)^2 + u(j+1)^2);
        endfor
    endfor

endfunction

# Algorytm optymalizacyjny
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, [],[],[],300);

obj
info
iter

t=[0:h:T]';

save wyniki.mat x l t N #zapis wyników

clear;
clc;


# Zmienne globalne
global number_of_parts = 3; #liczba członów
global number_of_time_steps=5; #liczba kroków czasowych
global T=1; #czas ruchu
global h=T/(number_of_time_steps-1); #krok całkowy
global mass=[1,1,1]; #masa
global partslength=[1,1.2,0.5]; #długość


# położenie końcowe końcówki układu
global xk=1;
global yk=1;


x0 = zeros(4*number_of_parts*number_of_time_steps,1);
# 4 ponieważ 4 wielkości angle, velocity, acceleration, torque


function equations = trapezoidal_transcription(x, xdot)

    global number_of_time_steps;
    global h;

    equations = zeros(number_of_time_steps-1,1);

    for i=[1:number_of_time_steps-1]
        equations(i) = x(i+1)-x(i)-1/2*h*(xdot(i+1)+xdot(i));
    endfor

endfunction


# Zdyskretyzowanie równania ruchu(równowagi)
function r = g (x)

    global number_of_time_steps;
    global number_of_parts;
    global h;
    global mass;
    global partslength;
    global xk;
    global yk;

    [angle, velocity, acceleration, torque] = unpacking_x(x, number_of_parts, number_of_time_steps);

    transcription_equations = []; boundary_conditions = [];

    for i=[1:number_of_parts]
        bufor = trapezoidal_transcription(angle(:,i), velocity(:,i));
        transcription_equations = [transcription_equations; bufor];
        bufor = trapezoidal_transcription(velocity(:,i), acceleration(:,i));
        transcription_equations = [transcription_equations; bufor];
        boundary_conditions = [boundary_conditions; angle(1,i); velocity(1,i); velocity(number_of_time_steps,i);];
    endfor

    boundary_conditions =[boundary_conditions; partslength*cos(angle(number_of_time_steps,:))'-xk; partslength*sin(angle(number_of_time_steps,:))'-yk]; # Warunki brzegowe
    
    motion_equations= creating_motion_equations(angle, velocity, acceleration, torque);

    r = [boundary_conditions; transcription_equations; motion_equations];

endfunction


# Funkcja celu - całka po czasie z sumy kwadratów momentów
function obj = phi (x)

    global number_of_time_steps;
    global number_of_parts;
    global h;
    obj = 0;

    for i=[1:number_of_parts]
        torque = x([3*number_of_time_steps*number_of_parts+(i-1)*number_of_time_steps+1:3*number_of_time_steps*number_of_parts+i*number_of_time_steps]);
        for j=[1:number_of_time_steps-1]
            obj = obj + 1/2 * h * (torque(j)^2 + torque(j+1)^2);
        endfor
    endfor

endfunction

# Algorytm optymalizacyjny
[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, [],[],[],300);

obj #końcowa wartość funkcji celu
info #jak się zakończył algorytm sqp: 101 - normalnie, 102 - error, 103 - wszystkie iteracje, 104 - krok stał się za mały
iter #liczba wykonanych iteracji przez sqp

t=[0:h:T]';

save wyniki.mat x partslength t number_of_parts #zapis wyników

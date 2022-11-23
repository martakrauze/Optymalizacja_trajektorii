clear;
clc;

load wyniki.mat

number_of_time_steps = rows(x)/(4*number_of_parts);

[p,v,a,u] = unpacking_x(x, number_of_parts, number_of_time_steps);

trajx = zeros(1,number_of_time_steps);
trajy = zeros(1,number_of_time_steps);

for i=[1:number_of_parts]
    trajx = [trajx; trajx(i,:) + (partslength(i) * cos(p(:,i)))'];
    trajy = [trajy; trajy(i,:) + (partslength(i) * sin(p(:,i)))'];
endfor

xkon=partslength*cos(p)';
ykon=partslength*sin(p)';

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
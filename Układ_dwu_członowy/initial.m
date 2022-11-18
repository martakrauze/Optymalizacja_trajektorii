clear;
clc;

# krok całkowania
global n=5;
global T=1;
global h=T/(n-1);

# przybliżenie początkowe
p10 = zeros(n,1);
v10 = zeros(n,1);
a10 = zeros(n,1);
u10 = zeros(n,1);
p20 = zeros(n,1);
v20 = zeros(n,1);
a20 = zeros(n,1);
u20 = zeros(n,1);

x = [p10; p20; v10; v20; a10; a20; u10; u20];

save wyniki.mat x
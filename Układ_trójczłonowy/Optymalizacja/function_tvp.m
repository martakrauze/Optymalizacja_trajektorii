function [x,xdot,xdotdot]=function_tvp(T,n,xp,xk)
    h=T/(n-1);
    t=[0:h:T];

    xdotc=1.5*(xk-xp)/T;
    xdotdotc=xdotc.^2./(xp-xk+xdotc*T);
    tc=(xp-xk+xdotc*T)./xdotc;
    x=zeros(n,1);
    xdot=zeros(n,1);
    xdotdot=zeros(n,1);

    for i=[1:n]
        if t(i)<=tc
        x(i)=xp+1/2*xdotdotc*t(i)^2;
        xdot(i)=xdotdotc*t(i);
        xdotdot(i)=xdotdotc;
        endif
        if t(i)>tc && t(i)<T-tc
        x(i)=xp+xdotdotc*tc*(t(i)-tc/2);
        xdot(i)=xdotc;
        xdotdot(i)=0;
        endif
        if t(i)>=T-tc
        x(i)=xk-1/2*xdotdotc*(T-t(i))^2;
        xdot(i)=xdotdotc*(T-t(i));
        xdotdot(i)=-xdotdotc;
        endif
    endfor
    
endfunction
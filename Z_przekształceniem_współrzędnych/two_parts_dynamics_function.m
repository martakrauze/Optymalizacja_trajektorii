function r=two_parts_dynamics_function(x,u)
    global I
    global Q 
    global M
    global grav
    global m1
    global m2
    global S1
    global S2

    u=[0;0;0];

   g=[0;grav*m1;u(1)-u(2);0;grav*m2;u(2)];
 
   thetadot=[x(3);x(4)];

   fi=[x(1);x(2)+x(1)];
   fidot=[x(3);x(3)+x(4)];

   S1_g=zeros(2,2);
   S2_g=zeros(2,2);
   for part=[1:1:2]
      S1_g(:,part)=R(fi(part))*S1(:,part);
      S2_g(:,part)=R(fi(part))*S2(:,part);
   endfor

   B=[d(1,1,S1_g,S2_g), zeros(2,1);
        1, 0;
       d(2,1,S1_g,S2_g), d(2,2,S1_g,S2_g);
          1, 1];  

   Bdot=[d_dot(1,1,S1_g,S2_g,fidot), zeros(2,1);
        0, 0;
       d_dot(2,1,S1_g,S2_g,fidot), d_dot(2,2,S1_g,S2_g,fidot);
          0, 0];

    M_new=B' * M * B;
    g_new=B' * (g - M * Bdot * thetadot);

    r1=x(3);
    r2=x(4);
    r3=M_new\g_new;
    r=[r1;r2;r3];

endfunction

function r=d(part, joint, S1_g, S2_g)
   global Q
   r=[0;0];
   for i=[joint:1:part]
      r=r-S1_g(:,i)+S2_g(:,i);
   endfor
   r=r-S2_g(:,part);
   r=Q*r;
endfunction

function r=d_dot(part, joint, S1_g, S2_g, fidot)
   global Q
   r=[0;0];
   for i=[joint:1:part]
      r=r+S1_g(:,i)*fidot(i)-S2_g(:,i)*fidot(i);
   endfor
   r=r+S2_g(:,part)*fidot(part);
endfunction

function rot=R(fi)
    rot=[cos(fi), -sin(fi); sin(fi), cos(fi)];
endfunction
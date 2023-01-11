function r=three_parts_dynamics_function_old(x,u)
    global I
    global Q 
    global M
    global s_A1
    global s_B1
    global s_B2
    global s_C2
    global s_C3
    global s_D3
    global grav
    global m1
    global m2
    global m3

    u=[0,0,0];
    g=[0;grav*m1;u(1)-u(2);0;grav*m2;u(2)-u(3);0;grav*m3;u(3)];

    fi1=x(1);
    fi2=x(2)+x(1);
    fi3=x(3)+x(2)+x(1);
    fidot1=x(4);
    fidot2=x(5)+x(4);
    fidot3=x(6)+x(5)+x(4);  
    thetadot=[x(4);x(5);x(6)];

    s_A1_g=R(fi1)*s_A1;
    s_B1_g=R(fi1)*s_B1;
    s_B2_g=R(fi2)*s_B2;
    s_C2_g=R(fi2)*s_C2;
    s_C3_g=R(fi3)*s_C3;
    s_D3_g=R(fi3)*s_D3;

    D=[I, Q * s_A1_g, zeros(2,3), zeros(2,3);
       I, Q * s_B1_g, -I, -Q * s_B2_g, zeros(2,3);
       zeros(2,3), I, Q * s_C2_g, -I, -Q * s_C3_g];
    

    d_11=-Q * s_A1_g;
    d_21=Q*(-s_A1_g+s_B1_g-s_B2_g);
    d_22=-Q*s_B2_g;
    d_31=Q*(-s_A1_g+s_B1_g-s_B2_g+s_C2_g-s_C3_g);
    d_32=Q*(-s_B2_g+s_C2_g-s_C3_g);
    d_33=-Q*s_C3_g;

    B=[d_11, zeros(2,1), zeros(2,1);
          1,          0,          0;
       d_21,       d_22, zeros(2,1);
          1,          1,          0;
       d_31,       d_32,       d_33;
          1,          1,         1];
    
    d_11_dot=s_A1_g*fidot1;
    d_21_dot=-(-s_A1_g*fidot1+s_B1_g*fidot1-s_B2_g*fidot2);
    d_22_dot=s_B2_g*fidot2;
    d_31_dot=-(-s_A1_g*fidot1+s_B1_g*fidot1-s_B2_g*fidot2+s_C2_g*fidot2-s_C3_g*fidot3);
    d_32_dot=-(-s_B2_g*fidot2+s_C2_g*fidot2-s_C3_g*fidot3);
    d_33_dot=s_C3_g*fidot3;

    Bdot=[d_11_dot, zeros(2,1), zeros(2,1);
                 0,          0,          0;
          d_21_dot,   d_22_dot, zeros(2,1);
                 0,          0,          0;
          d_31_dot,   d_32_dot,   d_33_dot,
                 0,          0,          0];
        

    M_new=B' * M * B;
    g_new=B' * (g - M * Bdot * thetadot);

    r1=x(4);
    r2=x(5);
    r3=x(6);
    r4=M_new\g_new;
    r=[r1;r2;r3;r4];

endfunction

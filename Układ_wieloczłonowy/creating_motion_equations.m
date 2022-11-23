function motion_equations= creating_motion_equations(angle, velocity, acceleration, torque)
    global number_of_time_steps;
    global number_of_parts;
    global mass;
    global partslength;

    motion_equations=[];

    for j=[1:number_of_time_steps]
        for i=[number_of_parts:-1:1]
            if i<number_of_parts
            joint_reaction_force(:,i)=fictitious_force(:,i+1)+joint_reaction_force(:,i+1);
            else
            joint_reaction_force(:,i)=zeros(2,1);
            endif
            fictitious_force(1,i)=-partslength(i)/2 * (cos(angle(j,i)) .* velocity(j,i).^2 + sin(angle(j,i)) .* acceleration(j,i));
            fictitious_force(2,i)=partslength(i)/2 * (-sin(angle(j,i)) .* velocity(j,i).^2 + cos(angle(j,i)) .* acceleration(j,i));
            vector=[-partslength(i)*sin(angle(j,i)), partslength(i)*cos(angle(j,i))];
            bufor=1/3 * mass(i) * partslength(i)^2 * acceleration(j,i) - torque(j,i) + vector * joint_reaction_force(:,i);
            motion_equations=[motion_equations;bufor];
        endfor
    endfor
    
endfunction
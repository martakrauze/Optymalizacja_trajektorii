function [angle, velocity, acceleration, torque]= unpacking_x(x, number_of_parts, number_of_time_steps)

    angle=[]; velocity=[]; acceleration=[]; torque=[];

    for i=[1:number_of_parts]
        bufor = x([(i-1)*number_of_time_steps+1:i*number_of_time_steps]);
        angle = [angle, bufor];
        bufor = x([number_of_time_steps*number_of_parts+(i-1)*number_of_time_steps+1:number_of_time_steps*number_of_parts+i*number_of_time_steps]);
        velocity = [velocity, bufor];
        bufor = x([2*number_of_time_steps*number_of_parts+(i-1)*number_of_time_steps+1:2*number_of_time_steps*number_of_parts+i*number_of_time_steps]);
        acceleration=[acceleration,bufor];
        bufor = x([3*number_of_time_steps*number_of_parts+(i-1)*number_of_time_steps+1:3*number_of_time_steps*number_of_parts+i*number_of_time_steps]);
        torque=[torque,bufor];
    endfor

endfunction
function dydt = ODE_function_no_lin_with_gravity(t,g,l,l0,k,m,y)
    %The non linearized EOM with gravity
    theta1 = y(1);
    theta2 = y(2);
    

    dydt = zeros(4,1);

    dydt(1) = y(3);
    dydt(2) = y(4);
    
    dydt(3:4) = Forward_dynamics_no_lin_with_gravity(g,k,l,l0,m,theta1,theta2);


end
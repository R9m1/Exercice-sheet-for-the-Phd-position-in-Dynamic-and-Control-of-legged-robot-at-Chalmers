function dydt = ODE_function(t,k,m,y)
    %The linearized EOM
    theta1 = y(1);
    theta2 = y(2);
    

    dydt = zeros(4,1);

    dydt(1) = y(3);
    dydt(2) = y(4);
    
    dydt(3:4) = Forward_dynamics(k,m,theta1, theta2);


end
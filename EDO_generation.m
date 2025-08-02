close all
clc

syms m l x l0 theta1 theta2 k dtheta1 dtheta2 theta10 theta20 g
syms t theta1_t(t) theta2_t(t) ddtheta1 ddtheta2 dtheta10 dtheta20

%Generation of the non linear EOM using Matlab symbolic toolbox
n = 1; %n=1 -> with gravity and n=2 ->without gravity
K1 = 1/2*m*l^2*dtheta1^2;
K2 = 1/2*m*l^2*dtheta2^2;
LT = sqrt((l0 + l*(sin(theta2) - sin(theta1)))^2 + (l*(cos(theta1)-cos(theta2)))^2);
Uspring = 1/2*k*(LT-l0)^2;
switch n
    case 1
        %The non linear EOM with gravity
        Ugrav = m*g*l*(cos(theta1) + cos(theta2));
    case 2
         %The non linear EOM without gravity
        Ugrav = 0;
end


Lag = K1 + K2 -Uspring -Ugrav

%First Lagrange equation of motion (theta1)
dLdtheta1_dot = diff(Lag,dtheta1);
dLdtheta1_dot = subs(dLdtheta1_dot, [dtheta1,dtheta2, theta1, theta2], [diff(theta1_t,t), diff(theta2_t,t),...
    theta1_t, theta2_t]);


dLdq1_dotdt = diff(dLdtheta1_dot,t);
dLdq1 = simplify(diff(Lag,theta1));
eq1 = dLdq1_dotdt - dLdq1;
eq1 = subs(eq1, [diff(theta1_t,t,t) diff(theta2_t,t,t) diff(theta1_t, t) diff(theta2_t, t) theta1_t theta2_t], [ddtheta1 ddtheta2 dtheta1 theta2 theta1 theta2]);


%Second Lagrange equation of motion (theta2)
dLdtheta2_dot = diff(Lag,dtheta2);
dLdtheta2_dot = subs(dLdtheta2_dot, [dtheta1,dtheta2, theta1, theta2], [diff(theta1_t,t), diff(theta2_t,t),...
    theta1_t, theta2_t]);

dLdq2_dotdt = diff(dLdtheta2_dot,t);
dLdq2 = simplify(diff(Lag,theta2));
eq2 = dLdq2_dotdt - dLdq2;
eq2 = subs(eq2, [diff(theta1_t,t,t) diff(theta2_t,t,t) diff(theta1_t, t) diff(theta2_t, t) theta1_t theta2_t], [ddtheta1 ddtheta2 dtheta1 theta2 theta1 theta2]);


eq =[eq1; eq2]
ddq = [ddtheta1; ddtheta2];
[M, C] = ODE_to_matrix_no_lin(eq,ddq);
qdd = M\(-C);


if n==1
    %With gravity
    matlabFunction(qdd, 'File', 'Forward_dynamics_no_lin_with_gravity', 'Optimize', false);
elseif n==2
    %Without gravity
     matlabFunction(qdd, 'File', 'Forward_dynamics_no_lin', 'Optimize', false);    
end 


%%
%Generation of the linear EOM
syms m l x l0 theta1 theta2 k dtheta1 dtheta2 theta10 theta20
syms t theta1_t(t) theta2_t(t) ddtheta1 ddtheta2 dtheta10 dtheta20


K1 = 1/2*m*l^2*dtheta1^2;
K2 = 1/2*m*l^2*dtheta2^2;
LT = l0 + l*(theta2  - theta1);
Uspring = 1/2*k*(LT-l0)^2;
Lag = K1 + K2 -Uspring;

%First Lagrange equation of motion (theta1)
dLdtheta1_dot = diff(Lag,dtheta1);
dLdtheta1_dot = subs(dLdtheta1_dot, [dtheta1,dtheta2, theta1, theta2], [diff(theta1_t,t), diff(theta2_t,t),...
    theta1_t, theta2_t]);


dLdq1_dotdt = diff(dLdtheta1_dot,t)
dLdq1 = simplify(diff(Lag,theta1));
eq1 = dLdq1_dotdt - dLdq1;
eq1 = subs(eq1, [diff(theta1_t,t,t) diff(theta2_t,t,t) diff(theta1_t, t) diff(theta2_t, t) theta1_t theta2_t], [ddtheta1 ddtheta2 dtheta1 theta2 theta1 theta2]);


%Second Lagrange equation of motion (theta2)
dLdtheta2_dot = diff(Lag,dtheta2);
dLdtheta2_dot = subs(dLdtheta2_dot, [dtheta1,dtheta2, theta1, theta2], [diff(theta1_t,t), diff(theta2_t,t),...
    theta1_t, theta2_t]);

dLdq2_dotdt = diff(dLdtheta2_dot,t);
dLdq2 = simplify(diff(Lag,theta2));
eq2 = dLdq2_dotdt - dLdq2;
eq2 = subs(eq2, [diff(theta1_t,t,t) diff(theta2_t,t,t) diff(theta1_t, t) diff(theta2_t, t) theta1_t theta2_t], [ddtheta1 ddtheta2 dtheta1 theta2 theta1 theta2]);


eq = simplify([eq1; eq2])

% eq = simplify([eq1; eq2])
ddq = [ddtheta1; ddtheta2];
q = [theta1; theta2];
q0 = [theta10; theta20];
[M,K,F] = ODE_to_matrix_lin(eq,ddq, q)
qdd = M\(-K*q)


matlabFunction(qdd,'File', 'Forward_dynamics', 'Optimize', false);
matlabFunction(M, 'File', 'Mass_matrix', 'Optimize',false);
matlabFunction(K, 'File', 'Stiffness_matrix', 'Optimize',false);
matlabFunction(F, 'File', 'F_Vector', 'Optimize',false);


qdd = subs(qdd, [theta1 theta2],[theta1_t theta2_t])

ode1 = diff(theta1_t,t,t) == qdd(1);
ode2 = diff(theta2_t,t,t) == qdd(2);
ode = [ode1; ode2]

dtheta1 = diff(theta1_t, t);
dtheta2 = diff(theta2_t, t);

conds = [ ...
    theta1_t(0) == 0;
    theta2_t(0) == -0.1;
    dtheta1(0) == 0;
    dtheta2(0) == 0;
];


q_ana = dsolve(ode,conds)
matlabFunction([q_ana.theta1_t; q_ana.theta2_t], 'File', 'Analytical_solution', 'Optimize',false);


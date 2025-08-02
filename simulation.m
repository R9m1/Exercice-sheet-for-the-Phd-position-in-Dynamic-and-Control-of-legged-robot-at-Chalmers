close all
clear all
%Task B, Problem 2 
l = 1;
l0 = 1;
m = 1;
k = 10;
tf = 5;
dt = 0.001;
y0 = [0;pi/3;0;0];
g = 9.81;


%Numerical simulation with ODE45 (4th order Runge-Kutta) for the linear EOM
tspan = 0:dt:tf;  
options = odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
disp("ODE45")
tic
[t,y] = ode45(@(t,y) ODE_function(t,k,m,y), tspan, y0,options);
close all
toc

theta1_RK45 = y(:,1);
theta2_RK45 = y(:,2);
theta1_dot_RK45 = y(:,3);
theta2_dot_RK45 = y(:,4); 



%Numerical simulation with the Newmark method for the linear EOM
nt = fix((tf)/dt);
M = Mass_matrix(l,m);
K = Stiffness_matrix(k,l);

F = zeros(length(M), 1);
C = zeros(length(M),length(M));

ti = 0;
beta = 1/4;
gamma = 1/6;
u0 = [y0(1); y0(2)];
du0 = [y0(3); y0(4)];
acceleration = 'Linear';
disp("Newmark")
tic
[dis, vel, acc, T] = Newmark_normalF(M,K,C,F,ti,tf,acceleration,u0,du0);
toc

theta1_newmark = dis(1,:);
theta2_newmark = dis(2,:);

%Analytical solution for the linear EOM
nt = fix((tf)/dt)
t = 0;
theta10 = y(1);
theta20 = y(2);
dtheta10 = y(3);
dtheta20 = y(4);

for i=1:nt+1
    sol = Analytical_solution(k,m,t);
    theta1_ana(i) = real(sol(1));
    theta2_ana(i) = real(sol(2));
    t = t + dt;
end


err_theta1_RK45 = abs(theta1_ana - theta1_RK45')./abs(theta1_ana);
err_theta2_RK45 = abs(theta2_ana - theta2_RK45')./abs(theta2_ana);
disp("Error RK45 for theta2")
disp(mean(err_theta2_RK45))

err_theta1_newmark = abs(theta1_ana - theta1_newmark)./abs(theta1_ana);
err_theta2_newmark = abs(theta2_ana - theta2_newmark)./abs(theta2_ana);
disp("Error Newmark for theta2")
disp(mean(err_theta2_newmark))


%Numerical simulation for the non linear EOM without gravity
tic
[t,y_no_lin] = ode45(@(t,y)  ODE_function_no_lin(t,l,l0,k,m,y), tspan, y0, options);
close all
toc

theta1_RK45_no_lin = y_no_lin(:,1);
theta2_RK45_no_lin = y_no_lin(:,2);


%Numerical simulation for the non linear EOM with gravity for three
%different initial conditions
tic
[t,y_no_lin_with_gravity1] = ode45(@(t,y)  ODE_function_no_lin_with_gravity(t,g,l,l0,k,m,y), tspan, y0, options);
close all
toc

y02 = [pi/2;pi/2;0;0];
y03 = [pi; pi; 0; 0];


[t,y_no_lin_with_gravity2] = ode45(@(t,y)  ODE_function_no_lin_with_gravity(t,g,l,l0,k,m,y), tspan, y02, options);
close all

[t,y_no_lin_with_gravity3] = ode45(@(t,y)  ODE_function_no_lin_with_gravity(t,g,l,l0,k,m,y), tspan, y03, options);
close all


theta1_RK45_no_lin_with_gravity1 = y_no_lin_with_gravity1(:,1);
theta2_RK45_no_lin_with_gravity1 = y_no_lin_with_gravity1(:,2);

theta1_RK45_no_lin_with_gravity2 = y_no_lin_with_gravity2(:,1);
theta2_RK45_no_lin_with_gravity2 = y_no_lin_with_gravity2(:,2);

theta1_RK45_no_lin_with_gravity3 = y_no_lin_with_gravity3(:,1);
theta2_RK45_no_lin_with_gravity3 = y_no_lin_with_gravity3(:,2);

%%
figure('Color', 'w')
subplot(1,2,1)
hold on
plot(tspan,theta1_RK45, "r", 'linewidth', 3); 
hold on
plot(T ,theta1_newmark, "b", 'linewidth', 3);
hold on
plot(tspan,theta1_ana, "g", 'linewidth', 3);
xlabel("Time t [s]")
ylabel("Angle \theta_1 [rad]")
legend("RK45","Newmark", "Analytical", 'Interpreter','latex',"Fontsize", 20)
grid on
set(gca, 'FontSize', 14)


subplot(1,2,2)
hold on
plot(tspan,theta2_RK45,  "r", 'linewidth', 3);
hold on
plot(T,theta2_newmark,  "b", 'linewidth', 3);
hold on
plot(tspan ,theta2_ana, "g", 'linewidth', 3)
xlabel("Time t [s]")
ylabel("Angle \theta_2 [rad]")
legend("RK45", "Newmark", "Analytical" , 'Interpreter','latex', "Fontsize", 20)
grid on
set(gca, 'FontSize', 14)


figure("Color", "w")
subplot(1,2,1)
hold on
plot(tspan,theta1_RK45, "r", 'linewidth', 3);
hold on
plot(tspan,theta1_RK45_no_lin, "b", 'linewidth', 3);
xlabel("Time t [s]")
ylabel("Angle \theta_1 [rad]")
legend("RK45 EOM linear", "RK45 EOM no linear", 'Interpreter','latex', "Fontsize", 15,...
    "location", "southwest")
grid on
set(gca, 'FontSize', 14)

subplot(1,2,2)
hold on
plot(tspan,theta2_RK45, "r", 'linewidth', 3);
hold on
plot(tspan,theta2_RK45_no_lin, "b", 'linewidth', 3);
xlabel("Time t [s]")
ylabel("Angle \theta_2 [rad]")
legend("RK45 EOM linear", "RK45 EOM non linear", 'Interpreter','latex',"Fontsize", 15,...
    "location", "southwest")
grid on
set(gca, 'FontSize', 14)


figure("Color", "w")
subplot(1,2,1)
semilogy(tspan,err_theta2_RK45, "r", 'linewidth', 3);
hold on
xlabel("Time t [s]")
ylabel("Error \theta_2 [rad]")
legend("Error RK45/analytical", 'Interpreter','latex', "Fontsize", 20, "Location", "southeast")
ylim([1e-11, 1e-3])
grid on
set(gca, 'FontSize', 14)

subplot(1,2,2)
semilogy(tspan,err_theta2_newmark, "k", 'linewidth', 3);
hold on
xlabel("Time t [s]")
ylabel("Error \theta_2 [rad]")
legend("Error Newmark/analytical", 'Interpreter','latex', "Fontsize", 20, "Location", "southeast")
ylim([1e-11, 1e-3])
grid on
set(gca, 'FontSize', 14)



% figure("Color","w")
% subplot(1,2,1)
% hold on
% plot(tspan, theta1_RK45_no_lin, "r", 'linewidth',3)
% hold on
% plot(tspan, theta1_RK45_no_lin_with_gravity, "k", 'linewidth',3)
% xlabel("Time t [s]");
% ylabel("Angle \theta_1 [rad]");
% legend("RK45 EOM non linear", "RK45 EOM non linear with gravity", 'Interpreter', 'latex', "Fontsize", 13, "location", "southwest");
% grid on;
% set(gca, 'FontSize', 14);
% 
% subplot(1,2,2)
% hold on
% plot(tspan, theta2_RK45_no_lin, "r", 'linewidth',3)
% hold on
% plot(tspan, theta2_RK45_no_lin_with_gravity, "k", 'linewidth',3)
% xlabel("Time t [s]");
% ylabel("Angle \theta_2 [rad]");
% legend("RK45 EOM non linear", "RK45 EOM non linear with gravity", 'Interpreter', 'latex', "Fontsize", 13, "location", "southwest");
% grid on;
% set(gca, 'FontSize', 14);


figure("Color","w")
subplot(1,2,1)
hold on
plot(tspan, theta1_RK45_no_lin_with_gravity1, "r", 'linewidth',3)
hold on
plot(tspan, theta1_RK45_no_lin_with_gravity2 , "k", 'linewidth',3)
hold on 
plot(tspan, theta1_RK45_no_lin_with_gravity3 , "b", 'linewidth',3)
hold on 
xlabel("Time t [s]");
ylabel("Angle \theta_1 [rad]");
legend("\theta_2(0) = \pi/3", "\theta_2(0) = \pi/2",  "\theta_2(0) = \pi",'Interpreter', 'latex', "Fontsize", 13, "location", "southwest");
grid on;
set(gca, 'FontSize', 14);

subplot(1,2,2)
hold on
plot(tspan, theta2_RK45_no_lin_with_gravity1, "r", 'linewidth',3)
hold on
plot(tspan, theta2_RK45_no_lin_with_gravity2 , "k", 'linewidth',3)
hold on 
plot(tspan, theta2_RK45_no_lin_with_gravity3 , "b", 'linewidth',3)
hold on 
xlabel("Time t [s]");
ylabel("Angle \theta_2 [rad]");
legend("\theta_2(0) = \pi/3", "\theta_2(0) = \pi/2",  "\theta_2(0) = \pi",'Interpreter', 'latex', "Fontsize", 13, "location", "southwest");
grid on;
set(gca, 'FontSize', 14);


%%
%Task B, Problem 3 (LQR controller)
A = [ 0, 0, 1, 0; 0, 0, 0, 1; ...
      -(g/l + k/m),  k/m, 0, 0; k/m, -(g/l + k/m), 0, 0];

B = [0; 0; 1/(m*l^2); 0];


Q = diag([1, 1, 1, 1]);  % penalize angles more than velocities
R = 0.01;                     % penalize torque

% Gain LQR
K_gain = lqr(A, B, Q, R);

x0 = [0; pi/3; 0; 0];  

% Simulation with ode45
% tspan = 0:dt:tf;
[t, x] = ode45(@(t, x) (A - B*K_gain)*x, tspan, x0);

% Plot
figure('Color','w');
plot(t, x,'linewidth',3);
xlabel('Temps (s)');
ylabel('States');
legend('\delta\theta_1', '\delta\theta_2', 'd\delta\theta_1', 'd\delta\theta_2');
% title('Response of the spring mechanism with the LQR control');
grid on;
set(gca, 'FontSize', 14);
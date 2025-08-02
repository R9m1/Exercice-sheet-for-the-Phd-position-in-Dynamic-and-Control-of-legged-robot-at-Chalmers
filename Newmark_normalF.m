function [dis,vel,acc,T] = Newmark_normalF(M,K,C,F,ti,tf,acceleration, u0, du0)
% NEWMARK'S METHOD : LINEAR SYSTEM
% Reference : Dynamics of Structures - Anil K. Chopra
%-------------------------------------------------------------------------
% Code written by :Siva Srinivas Kolukula                                 |
%                  Senior Research Fellow                                 |
%                  Structural Mechanics Laboratory                        |
%                  Indira Gandhi Center for Atomic Research               |
%                  INDIA                                                  |
% E-mail : allwayzitzme@gmail.com                                         |
%-------------------------------------------------------------------------
% Purpose : Dynamic Response of a system using linear Newmark's Method
% Synopsis :
%       [depl,vel,accl,U,t] = NewmarkMethod(M,K,C,P,phi,sdof,acceleration)
% 
% Variable Description :
% INPUT :
%       M - Mass Matrix (in modal coordinates)
%       K - Stiffness Matrix (in modal coordinates)
%       C - Damping Matrix (in modal coordinates)
%       P - Force Matrix (in modal coordinates)
%       sdof - system degree's of freedom
%       acceleration - Type of Newmark's Method to be used 
%       
% OUTPUT :
%        depl - modal displacement's 
%        vel - modal velocities
%        accl - modal accelerations 
%        U - system's displacement
%        t - time values at which integration is done
%--------------------------------------------------------------------------


switch acceleration
    case 'Average' 
        gaama = 1/2 ;beta = 1/4 ;
    case 'Linear'
        gaama = 1/2 ;beta = 1/6 ;
end
% Time step
dt = 0.001;
T = 0:dt:tf;
nt = fix((tf-ti)/dt);
n = length(M) ;

% Constants used in Newmark's integration

a0 = 1/(beta*dt^2) ;    a1 = gaama/(beta*dt) ;
a2 = 1/(beta*dt) ;       a3 = (1/(2*beta))-1 ;
a4 = (gaama/beta)-1 ;    a5 = dt/2*(gaama/beta-2) ;
a6=dt*(1-gaama);         a7 =dt*gaama;


dis = zeros(n,nt+1) ;
vel = zeros(n,nt+1) ;
acc = zeros(n,nt+1) ;
% Initial Conditions
dis(:,1) = u0 ;
vel(:,1) = du0 ;
P0=zeros(n,1) ;
acc(:,1) = inv(M)*(P0-C*vel(:,1)-K*dis(:,1)) ;



for i=1:nt
    
    ke=K+a0*M+a1*C;
    fe=F+M*(a0*dis(:,i)+a2*vel(:,i)+a3*acc(:,i))+C*(a1*dis(:,i)+a4*vel(:,i)+a5*acc(:,i));
    dis(:,i+1)=ke\fe;
    acc(:,i+1)=a0*(dis(:,i+1)-dis(:,i))-a2*vel(:,i)-a3*acc(:,i);
    vel(:,i+1)=vel(:,i)+a6*acc(:,i)+a7*acc(:,i+1);
    
    

end


end

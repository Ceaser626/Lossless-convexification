

clear all
clc


%------1¡¢Parameter Initialization
%---Background Parameter
g=[-3.7114 0 0];
m_dry=1505;
m_wet=1905;
I_sp=225;
T1=930;
T2=2480;
n=6;
phii=27/360*2*pi;
alpha=5.09e-4;
%---Optimal t_f and N
t_f=68;
N=55;
time_step=t_f/N;
%---Initial State
position_0=[3000 1000 2000];
velocity_0=[-50,10,100];
y0=[position_0 velocity_0 log(m_wet)]';%7*1


%------2¡¢Define Matrix
%---State Matrix
A_c=[zeros(3,3) eye(3) zeros(3,1);zeros(4,7)];
B_c=[zeros(3,3) zeros(3,1);eye(3) zeros(3,1);zeros(1,3) -alpha];
[A,B]=c2d(A_c,B_c,time_step);% A 7*7 ;B 7*4

%---Matrix for State Computation
phi=zeros(7,7*N); %A^k 7*7N
phi(:,1:7)=A;
for k=2:N
    phi(:,7*k-6:7*k)=A*phi(:,7*k-13:7*k-7);
end

lambda=zeros(7,4*N); %B+AB+......A^(k-1)B 7*4N
lambda(:,1:4)=B;
for k=2:N
    lambda(:,4*k-3:4*k)=A*lambda(:,4*k-7:4*k-4)+B;
end

psi=zeros(7*N,4*N+4); %psi 7N*(4N+4)
psi(1:7,1:4)=B;
for k=2:N
    psi(7*k-6:7*k,:)=A*psi(7*k-13:7*k-7,:);
    psi(7*k-6:7*k,4*k-3:4*k)=B;
end

upsilon=zeros(4*N,4*N+4); %upsilon 4N*(4N+4)
for k=1:N+1
    upsilon(4*k-3:4*k,4*k-3:4*k)=eye(4);
end

%---Matrix for Optimal Computation
Z0=log(m_wet-alpha*n*T2*cos(phii)*time_step*(0:N)'); %(N+1)*1
for k=0:N
    mu_1(k+1)=n*T1*cos(phii)*expm(-Z0(k+1));%(N+1)*1
    mu_2(k+1)=n*T2*cos(phii)*expm(-Z0(k+1));%(N+1)*1
end

epsilon_k=zeros(7,N);
for k=1:N
    epsilon_k(:,k)=phi(:,7*k-6:7*k)*y0+lambda(:,4*k-3:4*k)*[g 0]';
end

E_u=[eye(3) zeros(3,1)];
F=[zeros(1,6) 1];
e_sigma=[0 0 0 1];

omega=zeros(1,4*N); %1*4N
for k=1:N
    omega(1,4*k-3:4*k)=e_sigma;
end
omega=time_step*omega;


%------3¡¢Optimization
%---Solver
% cvx_solver  Mosek
cvx_solver  SDPT3

cvx_begin
%---Variable
variable eta(4*N+4,1);

%---Objective Function
minimize (omega*eta(1:4*N,1));

subject to
%---Final height constraint
[eye(3) zeros(3,4);zeros(4,3) zeros(4,4)]*(epsilon_k(:,N)+psi(7*N-6:7*N,:)*eta)==zeros(7,1);
%---Final velocity constraint
[0 0 0 1 0 0 0]*(epsilon_k(:,N)+psi(7*N-6:7*N,:)*eta)==0;
[0 0 0 0 1 0 0]*(epsilon_k(:,N)+psi(7*N-6:7*N,:)*eta)==0;
[0 0 0 0 0 1 0]*(epsilon_k(:,N)+psi(7*N-6:7*N,:)*eta)==0;

%---Convexified Thrust Constraint
for k=0:N
    norm(E_u*upsilon(4*k+1:4*k+4,:)*eta,2)<=e_sigma*upsilon(4*k+1:4*k+4,:)*eta;
end

%---thrust constraints
mu_1(1)*(1-(F*y0-Z0(1))+0.5*(F*y0-Z0(1)).^2)<=e_sigma*upsilon(1:4,:)*eta;
e_sigma*upsilon(1:4,:)*eta<=mu_2(1)*(1-(F*y0-Z0(1)));    
for k=1:N
    mu_1(k+1)*(1-(F*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta)-Z0(k+1))+0.5*(F*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta)-Z0(k+1)).^2)<=e_sigma*upsilon(4*k+1:4*k+4,:)*eta;
    e_sigma*upsilon(4*k+1:4*k+4,:)*eta<=mu_2(k+1)*(1-(F*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta)-Z0(k+1)));    
end

%---Fuel mass constraints
for k=1:N
    F*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta)<=log(m_wet-alpha*n*T1*cos(phii)*time_step*k);
    log(m_wet-alpha*n*T2*cos(phii)*time_step*k)<=F*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta);
end

% %---Cone constraints
% norm([0 1 0;0 0 1]*[eye(3) zeros(3,4)]*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta),2)+[-0.268 0 0 0 0 0 0]*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta)<=0;

cvx_end


%------Output Information
%Converting output into manageable format
display_position(1:3,N+1)=zeros;
display_velocity(1:3,N+1)=zeros;
display_netthrust(1:3,N+1)=zeros;
display_allthrust(1,N+1)=zeros;
display_allvelocity(1,N+1)=zeros;
display_position(1:3,1)=position_0;
display_velocity(1:3,1)=velocity_0;
mass(1)=m_wet;
for k=1:N
    display_position(1:3,k+1)=[eye(3) zeros(3,4)]*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta);
    display_velocity(1:3,k+1)=[zeros(3,3) eye(3) zeros(3,1)]*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta);
    mass(k+1)=[0 0 0 0 0 0 1]*(epsilon_k(:,k)+psi(7*k-6:7*k,:)*eta);
    mass(k+1)=exp(mass(k+1));
end
for k=0:N
    display_netthrust(1:3,k+1)=[eye(3) zeros(3,1)]*(upsilon(4*k+1:4*k+4,:)*eta);
    display_thrust(1:3,k+1)=display_netthrust(1:3,k+1);
    display_netthrust(1:3,k+1)=display_netthrust(1:3,k+1)*mass(k+1);
    display_allthrust(1,k+1)=norm(display_thrust(1:3,k+1),2)*mass(k+1);
    display_allvelocity(1,k+1)=norm(display_velocity(1:3,k+1),2);
end

%plot figures
figure
%3D Trajectory
subplot(3,3,7:9)
plot3(display_position(3, :), display_position(2, :), display_position(1, :), '-')
grid on
xlabel('z (meters)')
ylabel('y (meters)')
zlabel('x (meters)')
title('3D Trajectory')

%Position--Time
subplot(3,3,1)
plot(0:time_step:t_f,display_position)
grid on
legend('x','y','z');
xlabel('Time(s)')
ylabel('Position(m)')
title('Position--Time')

%Velocity--Time
subplot(3,3,2)
plot(0:time_step:t_f,display_velocity)
grid on
legend('x','y','z');
xlabel('Time(s)')
ylabel('Velocity(m/s)')
title('Velocity--Time')

%netthrust--Time
subplot(3,3,3)
plot(0:time_step:t_f,display_netthrust)
grid on
legend('x','y','z');
xlabel('Time(s)')
ylabel('Net thrust(N)')
title('netthrust--Time')

%Mass--Time
subplot(3,3,4)
plot(0:time_step:t_f,mass)
grid on
xlabel('Time(s)')
ylabel('Mass(kg)')
title('Mass--Time')

%allnetthrust--Time
subplot(3,3,6)
plot(0:time_step:t_f,display_allthrust)
grid on
xlabel('Time(s)')
ylabel('All Thrst(N)')
title('allnetthrust--Time')

%allvelocity--Time
subplot(3,3,5)
plot(0:time_step:t_f,display_allvelocity)
grid on
xlabel('Time(s)')
ylabel('All velocity(m/s)')
title('allvelocity--Time')

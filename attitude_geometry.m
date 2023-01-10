clear;clc;close all;
param;
delta_t = 0.01;
iteration = 1000; 


R = eye(3, 3);
R_dot = eye(3, 3);
e_R = zeros(3, 1);

omega = zeros(3, 1);
omega_d = zeros(3, 1);
e_Omega = zeros(3, 1);

J = diag([0.0820 0.0845 0.1377]);
alpha_record = zeros(iteration, 1);

R_init = angle2dcm(200*pi/180, 20*pi/180, 0); 
R = R_init;
R_d = eye(3, 3);
e_z_d = R_d(:, 3);

k_R = P.kR;
k_Omega = P.kOmega;


for i = 1:iteration
    e_z  = R(:, 3);
    alpha = acos(dot(e_z, e_z_d));
    alpha_record(i, :)= alpha;
    e_R = 1/2 * vee(R_d'*R - R'*R_d);
    e_Omega = omega - R'*R_d*omega_d;
    M = -k_R*e_R - k_Omega*e_Omega + cross(omega, J*omega);
    omega_dot = inv(J)*(M - cross(omega, J*omega));
    omega = omega + omega_dot*delta_t;
    
    R_dot = R*hat(omega);
    R = R + R_dot*delta_t;
    [U,S,V] = svd(R);
    R=U*V';
    R = real(R);
end

plot(1:iteration, alpha_record(1:iteration, 1));





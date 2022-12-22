clear;clc;close all;
delta_t = 0.01;
iteration = 1000;

% x:45 y:45 z:45
q_init = angle2quat(200*pi/180, 20*pi/180, 0); 

q = q_init;
q_record = zeros(1, 4,iteration);
q_dot = zeros(4, 1);
R = eye(3, 3);
R_dot = eye(3, 3);

eB_z = zeros(3, 1);
eB_z_des = [2 1 1];
eB_z_des = eB_z_des/norm(eB_z_des);

n = zeros(3, 1);
n_B = zeros(3, 1);
q_e_rp = zeros(1, 4);
alpha_degree = zeros(1, iteration);
gamma_degree = zeros(1, iteration);

omega_rp_des = zeros(1, 2);
omega_y_des = zeros(1, 1);
omega_des = zeros(1, 3);
omega = zeros(1, 3);
omega_dot = zeros(1, 3);

P = diag([40 40 10]);

yaw_des = 35*pi/180;

k_rp = 20;
k_y = 10;

J = diag([0.0820 0.0845 0.1377]);

for i = 1:iteration

    % Align roll, pitch
    eB_z = quatrotate(q, [0 0 1]);
    eB_x = quatrotate(q, [1 0 0]);
    eB_y = quatrotate(q, [0 1 0]);
    alpha = acos(dot(eB_z, eB_z_des));
    alpha_degree(1, i) = alpha*180/pi;
    n = cross(eB_z, eB_z_des)/norm(cross(eB_z, eB_z_des));
    n_B = quatrotate(quatinv(q), n);
    q_e_rp(1,:) = [cos(alpha/2) sin(alpha/2)*n_B];
    

    % Find omega_des for roll, pitch
    if q_e_rp(1) >= 0
        omega_rp_des = 2*k_rp*q_e_rp(2:3);
    else
        omega_rp_des = -2*k_rp*q_e_rp(2:3);
    end

    omega_des(1, :) = [omega_rp_des 0];

%     % Align yaw
%     if alpha*180/pi < 0.1
%         
%         eC_x = [cos(yaw_des); sin(yaw_des); 0];
%         eC_y = [-sin(yaw_des); cos(yaw_des); 0];
%         eB_x_des = cross(eC_y, eB_z_des)/norm(cross(eC_y, eB_z_des));
%         eB_y_des = cross(eB_z_des, eB_x_des)/norm(cross(eB_z_des, eB_x_des));
%         R = [eB_x_des' eB_y_des' eB_z_des'];
%         
%         q_des = rotm2quat(R);
%         q_e_y = quatmultiply(quatinv(quatmultiply(q, real(q_e_rp))), q_des);
%         gamma_degree(i) = acos(dot(eB_y, eB_y_des))*180/pi;
%         gamma_degree(i) = q_e_y(4);
%         
%         if q_e_y(1) >= 0
%             omega_y_des = 2*k_y*q_e_y(4);
%         else
%             omega_y_des = -2*k_y*q_e_y(4);
%         end
% 
%         omega_des(1, :) = [omega_rp_des omega_y_des];
%     end
    
    % Design controller
    tou = J*P*(omega_des(1, :)' - omega') + cross(omega', J*omega');
    omega_dot = inv(J)*(tou - cross(omega', J*omega'));
    omega = omega + omega_dot'*delta_t;
    
    % Calculate rotation
    R = quat2rotm(q);
    R_dot = R*hat(omega);
    R = R + R_dot*delta_t;
    [U,S,V] = svd(R);
    R=U*V';
    R = real(R);
    q = rotm2quat(R);
    
    q_record(1, :, i) = q;
    %e_R = 1/2 * vee(R_d'*R - R'*R_d);    

end
figure
plot(1:iteration, alpha_degree(1, :));
figure
plot(1:iteration, gamma_degree(1, :));
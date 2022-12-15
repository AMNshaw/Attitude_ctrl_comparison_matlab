function omega_rp_des = virtual_input(k, q_e_rp)

if q_e_rp(1) >= 0
    omega_rp_des = 2*k*q_e_rp(2:4);
else
    omega_rp_des = -2*k*q_e_rp(2:4);
end

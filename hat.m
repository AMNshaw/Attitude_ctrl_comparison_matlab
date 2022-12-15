function skew_matrix = hat(omega)

x = omega(1);
y = omega(2);
z = omega(3);

skew_matrix = [0 -z y;...
               z 0 -x;...
               -y x 0];
function XI = b_skew_se3(xi)
% compute the wedge ^ operator associated to se(3): xi^
% xi = [rho_x,rho_y,rho_z,phi_x,phi_y,phi_z]'

 XI = [b_skew_so3(xi(4:6)) xi(1:3);...
       0       0         0     0 ];    

end
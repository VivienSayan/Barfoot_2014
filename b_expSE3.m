function chi = b_expSE3(xi)
% xi (in R^6) [rho_x,rho_y,rho_z, rotation about x axis, rotation about y axis, rotation about z axis ]
% See Barfoot,2014 Appendix B - conversion from SE(3) to R^6
% See also Chirikjian Eq. 10.93
rho = xi(1:3);
phi = xi(4:6);
angles_norm = norm(phi);
if(angles_norm == 0)
    chi = [eye(3)        rho;...
           zeros(1,3)     1] ;
else
    Xi = [b_skew_so3(phi)       rho;... % in Lie algebra
           0     0     0         0];
    chi = eye(4) + Xi + 1/angles_norm^2*(1-cos(angles_norm))*Xi^2 + 1/angles_norm^3*(angles_norm-sin(angles_norm))*Xi^3; % in Lie group
end
end
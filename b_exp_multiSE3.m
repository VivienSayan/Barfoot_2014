function chi = b_exp_multiSE3(xi)
% xi (in R^{6+3*nbposes}) [x,y,z, rotation about x axis, rotation about y axis, rotation about z axis, x1,y1,z1, ... , xj, yj, zj]
% See "the invariant extended kalman filter as a stable observer,
% Barreau,Bonnabel, 2017" Appendix A.2
% See also https://www.ethaneade.com/lie.pdf (Eq.76-84)
rho = xi(1:3);
phi = xi(4:6);
add_rho = xi(7:end);
angles_norm = norm(phi);
Nb_rho = length([rho;add_rho])/3;
if(angles_norm == 0)
    chi = eye(3+Nb_rho);
    chi(1:3,4:end) = reshape([rho;add_rho],[3 Nb_rho]);
else
    Xi = zeros(3+Nb_rho);
    Xi(1:3,1:3) = b_skew_so3(rho);
    Xi(1:3,4:end) = reshape([rho;add_rho],[3 Nb_rho]);
    chi = eye(3+Nb_rho) + Xi + 1/angles_norm^2*(1-cos(angles_norm))*Xi^2 + 1/angles_norm^3*(angles_norm-sin(angles_norm))*Xi^3;
end

end
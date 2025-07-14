function Ad_chi = b_Adjoint_SE3(chi)
% See https://www.ethaneade.com/lie.pdf
R = chi(1:3,1:3);
T = b_skew_so3(chi(1:3,4));
Ad_chi = [R          T*R;...
          zeros(3,3)  R];
end
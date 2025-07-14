function inv_chi = b_invSE3(chi)
inv_chi = chi;
R = chi(1:3,1:3);
inv_R = R';
inv_chi(1:3,1:3) = inv_R;
inv_chi(1:3,4:end) = -inv_R*chi(1:3,4:end);
end


function [Rot,eul_angles,trans,add_pos] = b_chi2stateSE3(chi)
Rot = chi(1:3,1:3);
eul_angles = rotm2eul(Rot)'; % rotm2eul returns 1x3 array so transpose it to get 3x1
trans = chi(1:3,4);
if size(chi,2) > 4
    add_pos = chi(1:3,5:end);
else
    add_pos = [];
end
end

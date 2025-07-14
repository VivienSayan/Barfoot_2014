function ad_XI = b_ad_se3(XI)
% See Barfoot, 2014, (Eq.12) and (Eq.4)
rho = XI(1:3,4);
XI_phi = XI(1:3,1:3);
XI_rho = b_skew_so3(rho);
ad_XI = [XI_phi      XI_rho;...
          zeros(3,3) XI_phi];
end
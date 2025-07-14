clear all;
clc;

roll = 10*pi/180;
pitch = 20*pi/180;
yaw = 30*pi/180;
x1 = 0.1;
x2 = 0.2;
x3 = 0.3;
x = [x1;x2;x3];
eul_angles = [roll;pitch;yaw];
chi = b_state2chiSE3(eul_angles,x);

ksi_rotx = 1*pi/180;
ksi_roty = 2*pi/180;
ksi_rotz = 3*pi/180;
ksi_1 = 0.01;
ksi_2 = 0.02;
ksi_3 = 0.03;
ksi = [ksi_1;ksi_2;ksi_3;ksi_rotx;ksi_roty;ksi_rotz];
eta_rotx = -1*pi/180;
eta_roty = 4*pi/180;
eta_rotz = -3*pi/180;
eta_1 = -0.01;
eta_2 = 0.02;
eta_3 = -0.03;
eta = [eta_1;eta_2;eta_3;eta_rotx;eta_roty;eta_rotz];

% chi belongs to SE(3)
% ksi belongs to R^6
% unit test (Ad_chi * ksi)^ = chi * (ksi)^ * chi^(-1) 
% -> Barfoot Eq.108
Ad_chi = b_Adjoint_SE3(chi);
skew_Ad_chi_ksi = b_skew_se3(Ad_chi*ksi); % skew
res1 = skew_Ad_chi_ksi
skew_ksi = b_skew_se3(ksi);
res2 = chi*skew_ksi*b_invSE3(chi)

% unit test expm[(Ad_chi * ksi)^] = expSE3(Ad_chi*ksi) = chi * expm[(ksi)^] * chi^(-1) = chi * expSE3(ksi) * chi^(-1)
% -> Barfoot Eq.110
res11 = expm(skew_Ad_chi_ksi)
res11_ = b_expSE3(Ad_chi*ksi)
res22 = chi*expm(skew_ksi)*b_invSE3(chi)
res22_ = chi*b_expSE3(ksi)*b_invSE3(chi)

% unit test [KSI,ETA] = ( adse3(KSI)*eta )^
% -> Barfoot Eq.31
res111 = b_skew_se3(ksi)*b_skew_se3(eta)-b_skew_se3(eta)*b_skew_se3(ksi)
res222 = b_skew_se3( b_ad_se3(b_skew_se3(ksi)) * eta )

KSI = b_skew_se3(ksi);
Jl_SE3 = b_JacSE3(KSI,'RIGHT');
% properties |Jl_SE3| = |Jl_SO3|^2 -> Chirikjian vol.2 Eq. 10.96
res1111 = det(Jl_SE3)
res2222 = det(Jl_SE3(1:3,1:3))^2

% Chirikjian vol.2 Eq. 10.93
res4 = b_expSE3(ksi)
res5 = Jl_SE3(1:3,1:3)*ksi(1:3)
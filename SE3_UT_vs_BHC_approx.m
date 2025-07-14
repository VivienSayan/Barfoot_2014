clear; 
%close all; 
clc;

ERROR_SP = [];
ERROR_2nd = [];
ERROR_4th = [];
x_axis = linspace(0.01,1,20);
%rng(1);
for alpha = x_axis
alpha
ksi = [0;0;1;0;pi/4;0]; dim_state = length(ksi);
chi = b_expSE3(ksi);
P = alpha*diag([5; 10; 5; 1/2; 1; 1/2]);
S = chol(P,'lower');

omega = [0;2;0;pi/6;0;0];
OMEGA = b_expSE3(omega); 
Q = alpha*diag([10; 5; 5; 1/2; 1; 1/2]);
SQ = chol(Q,'lower');

chi_ = chi*OMEGA;
inv_chi_ = b_invSE3(chi_);

%------ UT -------------------------------------------------------------
% generer les sigma-points
P_aug = blkdiag(P,Q); dim_state_aug = size(P_aug,1); state_aug = zeros(dim_state_aug,1);
UT_kappa = 3-dim_state_aug;
[Wm,Wc,xi,XI,lambda] = b_compute_weights_SE3(dim_state_aug,1,0,UT_kappa); 
SigPts_01 = [zeros(dim_state_aug,1) -XI*eye(dim_state_aug,dim_state_aug) XI*eye(dim_state_aug,dim_state_aug)];
SigPts = zeros(dim_state_aug,2*dim_state_aug+1);
for i = 1:2*dim_state_aug+1
    SigPts(:,i) = state_aug(:) + sqrtm(P_aug)*SigPts_01(:,i);
end
SigPts_ut = SigPts; SigPts_qo = SigPts;
% Propagation
for j = 1:2*dim_state_aug+1
    ksi2_j = SigPts_ut(1:dim_state,j);
    ksi1_j = SigPts_ut(dim_state+1:end,j);
    chi2_j = chi * b_expSE3(ksi2_j);
    chi1_j = OMEGA * b_expSE3(ksi1_j);
    chi3_j = chi2_j * chi1_j;
    XI3_j = inv_chi_*chi3_j;
    SigPts_ut(1:dim_state,j) = b_logSE3(XI3_j);
end
%WSigPts_ = sqrt(Wc(1:end)).*SigPts_(1:dim_state,1:end);
%[~,RSx] = qr(WSigPts_');
%S_ = RSx(1:dim_state,1:dim_state);
%P_ = S_*S_'
P_ut = zeros(dim_state,dim_state);
for j = 1:2*dim_state_aug+1
    P_ut = P_ut + Wc(j)*SigPts_ut(1:dim_state,j)*SigPts_ut(1:dim_state,j)';
end


%----- P3 (2nd order) ----------------------------------------------------
U = b_invSE3( OMEGA );
Ad_SE3 = b_Adjoint_SE3(U);
AdPAd = Ad_SE3*P*Ad_SE3';
P_2nd = AdPAd + Q;

%----- P3 (4th order) ----------------------------------------------------
A1 = [linear_bracket(Q(4:6,4:6)) linear_bracket(Q(1:3,4:6)+Q(1:3,4:6)');...
       zeros(3,3)                linear_bracket(Q(4:6,4:6))];
A2p = [linear_bracket(AdPAd(4:6,4:6)) linear_bracket(AdPAd(1:3,4:6)+AdPAd(1:3,4:6)');...
       zeros(3,3)                     linear_bracket(AdPAd(4:6,4:6))];

Boo = bilinear_bracket(AdPAd(4:6,4:6),Q(4:6,4:6));

Bpo = bilinear_bracket(AdPAd(4:6,4:6),Q(1:3,4:6)') + bilinear_bracket(AdPAd(1:3,4:6)',Q(4:6,4:6));

Bpp = bilinear_bracket(AdPAd(4:6,4:6),Q(1:3,1:3)) + bilinear_bracket(AdPAd(1:3,4:6)',Q(1:3,4:6)) +...
      bilinear_bracket(AdPAd(1:3,4:6),Q(1:3,4:6)') + bilinear_bracket(AdPAd(1:3,1:3)',Q(4:6,4:6));

B = [Bpp  Bpo;...
     Bpo' Boo];
P_4th = P_2nd + 1/4*B + 1/12*(A2p*Q + Q*A2p' + A1*AdPAd + AdPAd*A1');

%----- Monte-Carlo -------------------------------------------------------
% generer les mc-points
M = 50000;
P_aug = blkdiag(P,Q); dim_state_aug = size(P_aug,1); state_aug = zeros(dim_state_aug,1);
MCPoints = sqrtm(P_aug)*randn(dim_state_aug,M);
% Propagation
for j = 1:M
    ksi2_j = MCPoints(1:dim_state,j);
    ksi1_j = MCPoints(dim_state+1:end,j);
    chi2_j = chi * b_expSE3(ksi2_j);
    chi1_j = OMEGA * b_expSE3(ksi1_j);
    chi3_j = chi2_j * chi1_j;
    XI3_j = inv_chi_*chi3_j;
    MCPoints(1:dim_state,j) = b_logSE3(XI3_j);
end
P_MC = zeros(dim_state,dim_state);
for j = 1:M
    P_MC = P_MC + 1/M*MCPoints(1:dim_state,j)*MCPoints(1:dim_state,j)';
end

error_sp = sqrt(trace((P_ut-P_MC)'*(P_ut-P_MC)));
error_2nd = sqrt(trace((P_2nd-P_MC)'*(P_2nd-P_MC)));
error_4th = sqrt(trace((P_4th-P_MC)'*(P_4th-P_MC)));

ERROR_SP = [ERROR_SP error_sp];
ERROR_2nd = [ERROR_2nd error_2nd];
ERROR_4th = [ERROR_4th error_4th];

end

save('tmp.mat')

figure()
plot(x_axis,ERROR_SP,'LineWidth',1); hold on; grid on;
plot(x_axis,ERROR_2nd,'--','LineWidth',1);
plot(x_axis,ERROR_4th,'--','LineWidth',1);
xlabel('$\alpha$','interpreter','latex');
ylabel('frobenius norm $P_{test}-P_{MC}$','interpreter','latex');
legend('$P_{sigma-point}$','$P_{2nd}$','$P_{4th}$','interpreter','latex');
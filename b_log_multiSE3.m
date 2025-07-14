function xi = b_log_multiSE3(chi)
% See Barfoot,2014 (Appendix B - conversion from SE(3) to R^6)
C = chi(1:3,1:3);
r = chi(1:3,4:end);
[vecs,eigs] = eig(C);
[~,idx] = min(abs(diag(eigs)-1));
a = vecs(:,idx); % instantaneous rotation axis

phi = acos(1/2*(trace(C)-1)); % instantaneous rotation angle; solution of tr(C) = 2cos(phi)+1
if phi == 0
    a = zeros(3,1);
    iJ = eye(3);
else
    A = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
    a_ta = a*transpose(a);
    % (Eq.97)
    C1 = cos(phi)*eye(3) + (1-cos(phi))*a_ta + sin(phi)*A;
    C2 = cos(-phi)*eye(3) + (1-cos(-phi))*a_ta + sin(-phi)*A;
    if norm(C2-C)<norm(C1-C)
        phi = -phi;
    end
    % (Eq.99)
    iJ = phi/2*cot(phi/2)*eye(3) + (1-phi/2*cot(phi/2))*a_ta - phi/2*A;
end
rho = iJ*r;
xi = real([rho(:);phi*a]);
end
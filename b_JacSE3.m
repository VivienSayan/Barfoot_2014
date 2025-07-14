function J_SE3 = b_JacSE3(XI,filter)
% See Barfoot,2014,"Associating uncertainty..." Appendix (Eq.100)
% RIGHT jacobian is for filter = LEFT equivariant-gaussian
% LEFT jacobian is for filter = RIGHT equivariant-gaussian

xi = [XI(1:3,4);XI(3,2);XI(1,3);XI(2,1)];

taille_xi = length(xi); % 3+3 + 3*Nb_add_pos
tolerance = 1e-1;
Nb_add_pos = taille_xi/3 -2; % =(1+1 + Nb_add_pos) -2
phi = xi(4:6);

PHI = XI(1:3,1:3);

ph = norm(phi);
if ph < tolerance % If the angle is small, fall back on the series representation
    switch filter
        case "LEFT"
            J_SO3 = b_JacSO3_Right_Series(PHI,10); % JacSO3 since phi is 3x1
        case "RIGHT"
            J_SO3 = b_JacSO3_Left_Series(PHI,10); % JacSO3 since phi is 3x1
        otherwise
            error("filter does not exist")
    end
    Q = b_xi2Q_Series(xi(1:6),10,filter);
    if Nb_add_pos > 0
        J_SE3 = kron(eye(Nb_add_pos+2),J_SO3);
        J_SE3(1:3,4:6) = Q;
        for i = 1:Nb_add_pos
            tmp = [xi(4+3*i:6+3*i); phi];
            Q = b_xi2Q_Series(tmp,10,filter);
            J_SE3(1:3,4+3*i:6+3*i) = Q; % :/ J_SE3(1+3*i:3+3*i,4+3*i:6+3*i)
        end
    else
        J_SE3 = [J_SO3          Q;...
                 zeros(3,3)  J_SO3];
    end

else % use closed form
    axis = phi/ph; % instantaneous rotation axis
    cph = (1 - cos(ph))/ph;
    sph = sin(ph)/ph;
    J_SO3 = sph * eye(3) + (1 - sph) * (axis)*(axis)' - cph * b_skew_so3(axis);
    Q = b_xi2Q(xi(1:6),filter);
    if Nb_add_pos > 0
        J_SE3 = kron(eye(Nb_add_pos+2),J_SO3);
        J_SE3(1:3,4:6) = Q;
        for i = 1:Nb_add_pos
            Q = b_xi2Q([xi(4+3*i:6+3*i);phi],filter);
            J_SE3(1:3,4+3*i:6+3*i) = Q; % :/ J_SE3(1+3*i:3+3*i,4+3*i:6+3*i)
        end
    else
        J_SE3 = [J_SO3          Q;...
                 zeros(3,3)  J_SO3];
    end

end

% must respect properties:
%det(J_SO3)^(2+Nb_add_pos) = det(J_SE3)

end
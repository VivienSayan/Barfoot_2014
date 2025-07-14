function J = b_JacSO3_Left_Series( XI, N )
% See Barfoot,2014 Appendix (Eq.98 in series representation)
xi = [XI(3,2);XI(1,3);XI(2,1)];
J = eye(3);
adn = eye(3);
ad = b_skew_so3(xi); % should be "little adjoint" operator but for so(3) we have adjoint(skew(xi)) = skew(xi)
for n = 1:N
    adn = adn*ad/(n + 1);    
    J = J + adn;
end

end
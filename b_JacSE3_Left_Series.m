function J = b_JacSE3_Left_Series( XI, N )
% See Barfoot,2014 Appendix (Eq.100 in series representation)
% See also LieAlg-UKF (Eq.4)
J = eye(6);
adn = eye(6);
ad = b_ad_se3(XI);
for n = 1:N
    adn = adn*ad /(n + 1);    
    J = J + adn;
end      

end
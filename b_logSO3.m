function res = b_logSO3(R)
% See https://www.ethaneade.com/lie.pdf (Eq.17-18)
theta = acos((trace(R)-1)/2) ;
if theta == 0
    tmp = zeros(3);
else
    tmp = theta/(2*sin(theta))*(R-R');
end
res = [-tmp(2,3);...
       tmp(1,3);...
       -tmp(1,2)];
end

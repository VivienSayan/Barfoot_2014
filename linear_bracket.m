function res = linear_bracket(A)
% Barfoot 2014 Eq.44

I = eye(size(A));
res = -trace(A)*I + A;

end
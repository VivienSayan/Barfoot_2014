function res = bilinear_bracket(A,B)
% Barfoot 2014 Eq.45

res = linear_bracket(A)*linear_bracket(B) + linear_bracket(A*B);

end
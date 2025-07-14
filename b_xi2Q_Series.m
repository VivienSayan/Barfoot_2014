function Q = b_xi2Q_Series(xi,N,filter)
% See Barfoot et al. (2014) Appendix (eq.102) 
switch filter
    case "RIGHT"
        rho = xi(1:3); % translation part
        phi = xi(4:6); % rotation part
    case "LEFT"
        rho = -xi(1:3); % translation part
        phi = -xi(4:6); % rotation part
    otherwise
        error("filter does not exist")
end

sum2 = 0;
for n=0:N
    sum1 = 1/factorial(n+2) * b_skew_so3(phi)^n * b_skew_so3(rho);
    for m=1:N
        sum1 = sum1 + 1/factorial(n+m+2) * b_skew_so3(phi)^n * b_skew_so3(rho) * b_skew_so3(phi)^m;
    end
    sum2 = sum2 + sum1;
end

Q = sum2;

end
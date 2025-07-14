function Q = b_xi2Q(xi,filter)
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

ph = norm(phi);
ph2 = ph*ph;
ph3 = ph2*ph;
ph4 = ph3*ph;
ph5 = ph4*ph;

cph = cos(ph);
sph = sin(ph);

rx = b_skew_so3(rho);
px = b_skew_so3(phi);

t1 = 0.5 * rx;
t2 = ((ph - sph)/ph3) * (px*rx + rx*px + px*rx*px);
m3 = ( 1 - 0.5 * ph2 - cph ) / ph4;
t3 = -m3 * (px*px*rx + rx*px*px - 3*px*rx*px);
m4 = 0.5 * ( m3 - 3*(ph - sph - ph3/6)/ph5 );
t4 = -m4 * (px*rx*px*px + px*px*rx*px);

Q = t1 + t2 + t3 + t4;

end
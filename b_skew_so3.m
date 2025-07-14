function XI = b_skew_so3(xi)
% compute the wedge ^ operator associated to so(3): xi^
% xi must be 3x1
XI = [0       -xi(3)   xi(2);...
      xi(3)     0     -xi(1);...
     -xi(2)   xi(1)     0];
end
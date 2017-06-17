function T = DH_transformation(alpha, d, theta, r)
% Modified DH params
T = [cos(theta)           -sin(theta)                     0              d;
     cos(alpha)*sin(theta) cos(alpha)*cos(theta) -sin(alpha) -r*sin(alpha);
     sin(alpha)*sin(theta) sin(alpha)*cos(theta)  cos(alpha)  r*cos(alpha);
                         0                     0           0            1];
end
%Small tolerance for special purposes
tol = 0.000001;
%Joint Angles
q1_lim = [-pi/2 pi/2];
q2_lim = [-pi/2 pi/2];
q3_lim = [-pi pi];
l2 = 0.20;
l3 = 0.15;

% The final point to be reached is given below
% The initial points are given in terms of joint coordinates.
qi = [0; 0; 0];
rfW = [0.05; 0.1; -0.1];
RotyPIby2 = [0 0 1;0 1 0;-1 0 0]; %Rotation from Robot to World Frame
rfR = inv(RotyPIby2)*rfW;
q = IGM_withConstraints(rfR(1), rfR(2), rfR(3), q1_lim, q2_lim, q3_lim, l2, l3)
Q = q*180/pi



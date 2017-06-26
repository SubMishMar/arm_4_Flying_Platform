function [rEF, rElbw] = DGM_func(q1, q2, q3, l2, l3)

% rEF is the position of end effector
% rElbw is the position of the Elbow
% Given below is wRr, the transform from robot frame to world frame.
RotyPIby2 = [0 0 1;0 1 0;-1 0 0];
T_RotyPIby2 = [RotyPIby2 zeros(3,1);0 0 0 1];

%DH Params
theta = [q1 q2 q3 0]; % Joint Angles
d = [ 0, 0, l2 l3]; 
r = [ 0, 0, 0, 0];
alpha = [ 0, pi/2, 0, 0];

T = zeros(4,4,5); %5 because we want to include the ground frame also
T(:,:,1) = eye(4,4);
T(:,:,2) = DH_transformation(alpha(1), d(1), theta(1), r(1));
for i = 3:5
     T(:,:,i) = T(:,:,i-1)*DH_transformation(alpha(i-1), d(i-1), theta(i-1), r(i-1));
end

% Given below is wRr, the transform from robot frame to world frame.
RotyPIby2 = [0 0 1;0 1 0;-1 0 0];


R_EF = RotyPIby2 * T(1:3,4,5);
R_Elbw = RotyPIby2 * T(1:3,4,4);

rEF = R_EF;
rElbw = R_Elbw;
end
clear 
clc

%% DH Parameters 
syms q1 q2 q3 q4
syms l2 l3 pi
theta = [q1 q2 q3 0]; % Joint Angles
d = [ 0, 0, l2 l3]; 
r = [ 0, 0, 0, 0];
alpha = [ 0, pi/2, 0, 0];

%% Transformation Matrices
T = sym(zeros(4,4,5));
T(:,:,1) = eye(4,4);
T(:,:,2) = simplify(DH_transformation(alpha(1), d(1), theta(1), r(1)));
for i = 3:5
     T(:,:,i) = simplify( T(:,:,i-1)*DH_transformation(alpha(i-1), d(i-1), theta(i-1), r(i-1)));
end


%% Jacobian Calculation

for i = 2:5
   az_i = T(1:3,1:3,i)*[0 ; 0 ; 1];
   w = T(1:3,4,5) - T(1:3,4,i);
   J(:,i-1) = [cross(az_i,w); az_i];
end
J = simplify(J);
J = J(1:6,1:3)
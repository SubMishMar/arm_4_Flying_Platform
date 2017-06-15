function [x, y, z] = DGM(q1, q2, q3, l2, l3)
% Z axis is along the first joint (or lets call it the pitching joint)
% (x, y, z) is the position of the end effector
x = ( l2*cos(q2) + l3*cos(q2+q3) )*cos(q1);
y = ( l2*cos(q2) + l3*cos(q2+q3) )*sin(q1);
z = ( l2*sin(q2) + l3*sin(q2+q3) );
end
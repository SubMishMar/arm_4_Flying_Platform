%clear all
%clc

l2 = 0.20;
l3 = 0.20;



xW = 0.1;
yW = 0.1;
zW = -0.1;

RotyPIby2 = [0 0 1;0 1 0;-1 0 0];
rR = inv(RotyPIby2)*[xW; yW; zW];
xR = rR(1);
yR = rR(2);
zR = rR(3);
[xR, yR, zR];
[q1, q2, q3] = IGM(xR, yR, zR, l2, l3);

[q1;q2;q3;]
Q1 = q1*180/pi;
Q2 = q2*180/pi;
Q3 = q3*180/pi;
[Q1;Q2;Q3;]

if isreal(q1) && isreal(q2) && isreal(q3) 
    
    plotter(q1(1), q2(1), q3(1), l2, l3, 1, 0)
    plotter(q1(2), q2(2), q3(2), l2, l3, 2, 0)
    plotter(q1(3), q2(3), q3(3), l2, l3, 3, 0)
    plotter(q1(4), q2(4), q3(4), l2, l3, 4, 0)

else
    
    disp(['Enter a Valid Point in the Workspace']);

end

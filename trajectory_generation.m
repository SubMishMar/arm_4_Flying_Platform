%Joint Angles
q1_lim = [-pi/2 pi/2];
q2_lim = [-pi/2 pi/2];
q3_lim = [-pi pi];
l2 = 0.20;
l3 = 0.40;

% The final point to be reached is given below
% The initial points are given in terms of joint coordinates.
qi = [0; 0; 0];
rfW = [-0.1; 0.2; -0.3];
RotyPIby2 = [0 0 1;0 1 0;-1 0 0];
rfR = inv(RotyPIby2)*rfW;
[qf1, qf2, qf3] = IGM(rfR(1), rfR(2), rfR(3), l2, l3)
qf_accept = [];
for i=1:4
    if qf1(i) >= q1_lim(1) && qf1(i) <= q1_lim(2)
       if qf2(i) >= q2_lim(1) && qf2(i) <= q2_lim(2)
           if qf3(i) >= q3_lim(1) && qf3(i) <= q3_lim(2)
               qf_accept = [qf_accept; qf1(i) qf2(i) qf3(i)];
           end
       end
    end    
end



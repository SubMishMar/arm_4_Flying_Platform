function qf_accept = IGM_withConstraints(x, y, z, q1_lim, q2_lim, q3_lim, l2, l3)
% qf_accept: accepted joint configs, each row is a set of joint angles that
% satisfies the constraints
%tol = 0;% tunable tolerance
[qf1, qf2, qf3] = IGM(x, y, z, l2, l3);
if (isreal(qf1) && isreal(qf2) && isreal(qf3))
qf_accept = [];
prev_elb = 0; %previous elbow z coordinates
l = 1;
for i = 1 : 4
    [qf1(i)*180/pi; qf2(i)*180/pi; qf3(i)*180/pi];
    if (qf1(i) >= (q1_lim(1))) && (qf1(i) <= (q1_lim(2)))
       if (qf2(i) >= (q2_lim(1))) && (qf2(i) <= (q2_lim(2)))
           if (qf3(i) >= (q3_lim(1))) && (qf3(i) <= (q3_lim(2)))
               [rEF, rElbw] = DGM_func(qf1(i), qf2(i), qf3(i), l2, l3);
               %[rElbw(3),prev_elb]
               %[qf1(i) qf2(i) qf3(i)]
               if rElbw(3) < 0 && rElbw(3)<=prev_elb  % Elbow shud be below the plane
                 if rElbw(3) < prev_elb
                    qf_accept(l,:) = [qf1(i) qf2(i) qf3(i)];
                 end
                 if rElbw(3) == prev_elb
                   qf_accept = [qf_accept; qf1(i) qf2(i) qf3(i)];
                   l = l+1;
                 end
               prev_elb = rElbw(3);
               end
           end
       end
    end    
end
nrows = size(qf_accept, 1);

if (~isempty(qf_accept))
  for m=1:nrows
      plotter(qf_accept(m,1), qf_accept(m,2), qf_accept(m,3), l2, l3, m, 0)
  end
else
  qf_accept = [];
  disp('No Suitable Solution [type 2]');
end

else
  qf_accept = [];
  disp('No Suitable Solution [type 1]');
end
end
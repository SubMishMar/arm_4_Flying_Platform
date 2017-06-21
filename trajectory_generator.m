function W = trajectory_generator(r_i, r_f, q1_lim, q2_lim, q3_lim, l2, l3, Tsim)
% Steps
% 1. Idea is to take initial and final task space coordinates as input, then
% 2. pass it to the IGM to find the corresponding joint space coordinates
% 3. find the way points in the joint space coordinates
% 4. Apply DGM to each of these way points to find the corresponding point in
%    the task space, use this task space coordinate via points for the
%    trajectory generator in the controller
[q1_i, q2_i, q3_i] = IGM(r_i(1), r_i(2), r_i(3), l2, l3);
[q1_f, q2_f, q3_f] = IGM(r_f(1), r_f(2), r_f(3), l2, l3);

end
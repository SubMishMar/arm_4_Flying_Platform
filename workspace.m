function W=workspace(q1_lim, q2_lim, q3_lim, l2, l3, res)
W=[];
 for q1 = linspace(q1_lim(1), q1_lim(2), res)
     for q2 = linspace(q2_lim(1), q2_lim(2), res)
         for q3 = linspace(q3_lim(1), q3_lim(2), res)
            [rEF, rElbw] = DGM_func(q1, q2, q3, l2, l3);
             W = [W;rEF'];
         end
     end
 end
end
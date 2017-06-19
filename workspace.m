function W=workspace(q1_lim, q2_lim, q3_lim, l2, l3)
W=[];
 for q1 = q1_lim(1):0.25:q1_lim(2)
     for q2 = q2_lim(1):0.25:q2_lim(2)
         for q3 = q3_lim(1):0.25:q3_lim(2)
            [x,y,z] = DGM(q1, q2, q3, l2, l3);
             W = [W;x y z];
         end
     end
 end
end
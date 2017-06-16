syms q1 q2 q3 l2 l3

J = [-sin(q1)*(l2*cos(q2)+l3*cos(q2+q3)) -cos(q1)*(l2*sin(q2)+l3*sin(q2+q3)) -l3*cos(q1)*sin(q2+q3);
      cos(q1)*(l2*cos(q2)+l3*cos(q2+q3)) -sin(q1)*(l2*sin(q2)+l3*sin(q2+q3)) -l3*sin(q1)*sin(q2+q3);
                                      0   l2*cos(q2)+l3*cos(q2+q3)            l3*cos(q2+q3)];
determinant = simplify(det(J))
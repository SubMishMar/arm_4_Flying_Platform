l2 = 0.20;
l3 = 0.15;

syms q2 q3 pi

ezplot(-l2*l3*(l2*cos(q2)*sin(q3) - l3*sin(q2) + l3*cos(q3)^2*sin(q2) + l3*cos(q2)*cos(q3)*sin(q3)) == 0)
title('Singular Points')



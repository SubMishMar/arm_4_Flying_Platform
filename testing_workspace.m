q1_lim = [-pi/2 pi/2];
q2_lim = [-pi/2 pi/2];
q3_lim = [-pi pi];
l2 = 0.15;
l3 = 0.15;
res = 20;
W=workspace(q1_lim , q2_lim, q3_lim, l2, l3, res);
plot3(W(:,1),W(:,2),W(:,3),'.')
axis square
grid


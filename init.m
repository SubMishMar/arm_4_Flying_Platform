%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MDH Table-- Refer to the book by Wisama Khalil
%on Modelling, Control and Identification of 
%Robots
%aplha    |  d   |theta |  r  |
%------------------------------ 
%  0      |  0   | q1   |  0  |
%  pi/2   |  0   | q2   |  0  |
%  0      |  l2  | q3   |  0  |
%  0      |  l3  | 0    |  0  |

% Joint limits
% q1_lim = [-pi/2 pi/2];
% q2_lim = [-pi/2 pi/2];
% q3_lim = [-pi pi];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation
l2 = 0.20;
l3 = 0.20;
K = 50*diag(ones(1,3));
q1_lim = [-pi/2 pi/2];
q2_lim = [-pi/2 pi/2];
q3_lim = [-pi pi];
x_i = 0.3; y_i = 0; z_i = -0.1;
x_f = 0.1; y_f = 0.1; z_f = -0.2 ;
tf = 100; %time of trajectory
sim('kinematic_model_control');

%% Plotter
nrows_q = size(q, 1);
figure(1);
i=1;
plotter(q(i,1), q(i,2), q(i,3), l2, l3, 1, 0);
hold on;
for i=1:100:nrows_q
plot3(r(:,1),r(:,2),r(:,3),'-mo','MarkerFaceColor',[1 0 0],'MarkerSize',1,'MarkerEdgeColor','r');
hold on;
plotter(q(i,1), q(i,2), q(i,3), l2, l3, 1, 0);
pause(0.01);
end

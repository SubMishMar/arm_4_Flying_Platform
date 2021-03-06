function plotter(q1, q2, q3, l2, l3, n, hold_on_off)
% n is the figure number
% hold_on_off to turn hold on = 1 or off = any other number
%% DH Params for the joint config
theta = [q1(1) q2(1) q3(1) 0];
r = zeros(1,4);
alpha = [ 0, pi/2, 0, 0];
d = [ 0, 0, l2 l3]; 


%% Transformation Matrices
T = zeros(4,4,5); %5 because we want to include the ground frame also
T(:,:,1) = eye(4,4);
T(:,:,2) = DH_transformation(alpha(1), d(1), theta(1), r(1));
for i = 3:5
     T(:,:,i) = T(:,:,i-1)*DH_transformation(alpha(i-1), d(i-1), theta(i-1), r(i-1));
end

% Given below is wRr, the transform from robot frame to world frame.
RotyPIby2 = [0 0 1;0 1 0;-1 0 0];
%RotyPIby2 = eye(3,3);
T_RotyPIby2 = [RotyPIby2 zeros(3,1);0 0 0 1];

%% Plotter
figure(n)
for i = 1:5
   a_xi = T_RotyPIby2*T(:,:,i)*[ 0.05 ; 0 ; 0 ; 1 ]; a_xi = a_xi(1:3);
   a_yi = T_RotyPIby2*T(:,:,i)*[ 0 ; 0.05 ; 0 ; 1 ]; a_yi = a_yi(1:3);
   a_zi = T_RotyPIby2*T(:,:,i)*[ 0 ; 0 ; 0.05 ; 1 ]; a_zi = a_zi(1:3);
   P_i = RotyPIby2*T(1:3,4,i);
   if i>1
       P_i_1 = RotyPIby2*T(1:3,4,i-1);
       % Plottting Links
       plot3(linspace(P_i_1(1),P_i(1),100),linspace(P_i_1(2),P_i(2),100),....
       linspace(P_i_1(3),P_i(3),100),'c', 'LineWidth',4);
       hold on;
   end
    % Transformed Frames
    if i==1
    plot3(linspace(P_i(1),a_xi(1),100),linspace(P_i(2),a_xi(2),100),....
          linspace(P_i(3),a_xi(3),100),'r','LineWidth',3)
    hold on;
    plot3(linspace(P_i(1),a_yi(1),100),linspace(P_i(2),a_yi(2),100),....
          linspace(P_i(3),a_yi(3),100),'g','LineWidth',3)
    hold on;
    plot3(linspace(P_i(1),a_zi(1),100),linspace(P_i(2),a_zi(2),100),....
          linspace(P_i(3),a_zi(3),100),'b','LineWidth',3)
    hold on;
    else
        plot3(linspace(P_i(1),a_xi(1),100),linspace(P_i(2),a_xi(2),100),....
          linspace(P_i(3),a_xi(3),100),'r','LineWidth',1)
    hold on;
    plot3(linspace(P_i(1),a_yi(1),100),linspace(P_i(2),a_yi(2),100),....
          linspace(P_i(3),a_yi(3),100),'g','LineWidth',1)
    hold on;
    plot3(linspace(P_i(1),a_zi(1),100),linspace(P_i(2),a_zi(2),100),....
          linspace(P_i(3),a_zi(3),100),'b','LineWidth',1)
    hold on;
    end
end
grid;
xlim([-1.2*(l3+l2) 1.2*(l3+l2)]);
ylim([-1.2*(l3+l2) 1.2*(l3+l2)]);
zlim([-1.2*(l3+l2) 1.2*(l3+l2)]);
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
title(['Simplified schematic of the aerial arm in joint config ', num2str(n)])
pbaspect([1 1 1])
axis square
if hold_on_off == 1
    
    hold on;
    
else
    
    hold off;
    
end

end
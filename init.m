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

l2 = 0.40;
l3 = 0.20;
K = 5*diag(ones(1,3));
sim('kinematic_model_control');
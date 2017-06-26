function [q1, q2, q3] = IGM(x, y, z, l2, l3)
% q1(1), q2, q3 store the corresponding joint angles 
% in radians
% x, y, z are in robot base frame
q1 = zeros(1, 4);
q2 = zeros(1, 4);
q3 = zeros(1, 4);

if x<0 && y == 0
  q1(1) = 0;
  q1(2) = q1(1);
  q1(3) = pi;
  q1(3) = atan2(sin(q1(3)), cos(q1(3)));
  q1(4) = q1(3);
elseif x==0 && y<0 
  q1(1) = pi/2;
  q1(2) = q1(1);
  q1(3) = 3*pi/2;
  q1(3) = atan2(sin(q1(3)), cos(q1(3)));
  q1(4) = q1(3);
else
  q1(1) = atan2(y,x);
  q1(2) = q1(1);
  q1(3) = atan2(y,x) + pi;
  q1(3) = atan2(sin(q1(3)), cos(q1(3)));
  q1(4) = q1(3);   
end

if sin(q1(1)) >= -0.000001 && sin(q1(1)) <= 0.000001
    
    %disp('sin(q1(1)) is zero');

    q3(1) = acos((x^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    if(isreal(q3(1))) %% not a smart way to check, can also check if the
        %point is within the workspace by checking if the thing within the 
        % acos() function is between -1 and 1, but whatever.. 

    q2(1) = (atan2(z, x) - atan2(l3*sin(q3(1)), l2 + l3*cos(q3(1))));
 

    else
        disp('Point Not Reachable');
        return;
    end
  
    q3(2) = -acos((x^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    if(isreal(q3(2)))

    q2(2) = atan2(z, x) - atan2(l3*sin(q3(2)), l2 + l3*cos(q3(2)));   

    else
        disp('Point Not Reachable');
        return;
    end
    
elseif cos(q1(1)) >= -0.000001 && cos(q1(1)) <= 0.000001
    
    %disp('cos(q1(1)) is zero');
    
    q3(1) = acos((y^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    
    if(isreal(q3(1)))

    q2(1) = atan2(z, y) - atan2(l3*sin(q3(1)), l2 + l3*cos(q3(1))) ;

    else
        disp('Point Not Reachable');
        return;
    end
    
    q3(2) = -acos((y^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    if(isreal(q3(2)))

    q2(2) = atan2(z, y) - atan2(l3*sin(q3(2)), l2 + l3*cos(q3(2)));   

    else
        disp('Point Not Reachable');
        return;
    end
    
else
    %disp('Neither sin(q1(1)) nor cos(q1(1)) is zero');
    
    xhat = x/cos(q1(1));
    
    q3(1) = acos((xhat^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    if(isreal(q3(1)))
    q2(1) = atan2(z, xhat) - atan2(l3*sin(q3(1)), l2 + l3*cos(q3(1)));
    else
        disp('Point Not Reachable');
        return;
    end
    
    q3(2) = -acos((xhat^2 + z^2 - l2^2 -l3^2)/(2*l2*l3)); 
    if(isreal(q3(2)))
    q2(2) = atan2(z, xhat) - atan2(l3*sin(q3(2)), l2 + l3*cos(q3(2)));
    else
       disp('Point Not Reachable');
       return;
    end
end


q2(1) = atan2(sin(q2(1)), cos(q2(1)));
q2(2) = atan2(sin(q2(2)), cos(q2(2)));
q2(3) = (pi-q2(1));
q2(3) = atan2(sin(q2(3)), cos(q2(3)));
q2(4) = (pi-q2(2));
q2(4) = atan2(sin(q2(4)), cos(q2(4)));

q3(1) = atan2(sin(q3(1)), cos(q3(1)));
q3(2) = atan2(sin(q3(2)), cos(q3(2)));
q3(3) = -q3(1);
q3(4) = -q3(2);


end
function [q1, q2, q3] = IGM(x, y, z, l2, l3)
q1 = 0;
q2 = zeros(1, 2);
q3 = zeros(1, 2);


q1 = atan2(y, x);

    
if sin(q1) >= -0.000001 && sin(q1) <= 0.000001
    
    %disp('sin(q1) is zero');
    
    q3(1) = acos((x^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    if(isreal(q3(1))) %% not a smart way to check, can also check if the
        %point is within the workspace by checking if the thing within the 
        % acos() function is between -1 and 1, but whatever.. 
    if x < 0
    q2(1) = atan2(z, x) - atan2(l3*sin(q3(1)), l2 + l3*cos(q3(1))) + pi;
    else
    q2(1) = atan2(z, x) - atan2(l3*sin(q3(1)), l2 + l3*cos(q3(1))) ;
    end
    else
        disp('Point Not Reachable');
        return;
    end
    
    q3(2) = -acos((x^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    if(isreal(q3(2)))
    if x < 0
    q2(2) = atan2(z, x) - atan2(l3*sin(q3(2)), l2 + l3*cos(q3(2))) + pi;
    else
    q2(2) = atan2(z, x) - atan2(l3*sin(q3(2)), l2 + l3*cos(q3(2)));   
    end
    else
        disp('Point Not Reachable');
        return;
    end
    
elseif cos(q1) >= -0.000001 && cos(q1) <= 0.000001
    
    %disp('cos(q1) is zero');
    
    q3(1) = acos((y^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    
    if(isreal(q3(1)))
    if y < 0
    q2(1) = atan2(z, y) - atan2(l3*sin(q3(1)), l2 + l3*cos(q3(1)))+pi;
    else
    q2(1) = atan2(z, y) - atan2(l3*sin(q3(1)), l2 + l3*cos(q3(1))) ;
    end
    else
        disp('Point Not Reachable');
        return;
    end
    
    if(isreal(q3(2)))
    q3(2) = -acos((y^2 + z^2 - l2^2 -l3^2)/(2*l2*l3));
    if y < 0
    q2(2) = atan2(z, y) - atan2(l3*sin(q3(2)), l2 + l3*cos(q3(2)))+pi;
    else
    q2(2) = atan2(z, y) - atan2(l3*sin(q3(2)), l2 + l3*cos(q3(2)));   
    end
    else
        disp('Point Not Reachable');
        return;
    end
    
else
    %disp('Neither sin(q1) nor cos(q1) is zero');
    
    xhat = x/cos(q1);
    
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
end
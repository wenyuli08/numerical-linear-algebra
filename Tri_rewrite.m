% function for 7.3.6
function [R] = Tri_rewrite(b)
    d(1:50) = 5;
    if length(b) ~= 50
        disp("A is a 50x50 tridiagonal matrix, so b must be (50,1).")
    end
    
    for i = 2:50
        % xmult = -1/5
        d(i) = d(i) + 1/5;
        b(i) = b(i) + (1/5)*b(i-1);
    end
    
    % Back substitution
    R(50) = b(50)/d(50);
    for i = 49:-1:1
        R(i) = (b(i) - R(i+1))/d(i);
    end
end
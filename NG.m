function [R_NG] = NG(A2,b2)
    [m, n] = size(A2);
    [rb, cb] = size(b2);
    R_NG = zeros(m,1);
    
    % data preparation
    if m ~= n
        disp("A must be a square matrix.");
    end
    if rb ~= m || cb ~= 1
        disp("b must be a m x 1 matrix.");
    end
    d = det(A2);
    if d == 0
        disp("Warning: singular matrices have either none or infinite number of solutions.")
    end
    
    % Forward elimination
    for k = 1:(n-1)
        for i = (k+1):m
            % the first(k-th) row is the pivot, ak,k is the pivot element
            xmult = A2(i,k)/A2(k,k);
            A2(i,k) = xmult;
            % calculate the (k+1)-th row * xmult
            for j = (k+1):n
                A2(i,j) = A2(i,j) - xmult*A2(k,j);
            end
            b2(i) = b2(i) - xmult*b2(k);
        end
    end
    
    disp(A2)
    disp(b2)
    
    % Back substitution
    R_NG(m) = b2(m)/A2(m,n);
    for i = (n-1):-1:1
        sum = b2(i);
        for j = (i+1):n
            sum = sum - A2(i,j)*R_NG(j);
        end
        R_NG(i) = sum/A2(i,i);
    end
end
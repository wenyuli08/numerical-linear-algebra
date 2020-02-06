function [R] = Hessenberg(A,b)
    [m,n] = size(A);
    [rb, cb] = size(b);
    if (m ~= n)
        disp("Hessenberg matrix should be a square matrix.")
    end
    if rb ~= m || cb ~= 1
        disp("Invalid b")
    end
    disp([A,b])
    
    for i = 1:(m-1)
        % deal with 0s on the main diagonal
        if A(i,i) == 0
            A(i,i) = 0.0001;
        end
        % perform Naive Gaussian
        xmult = A(i+1,i)/A(i,i);
        A(i+1,:) = A(i+1,:) - xmult*A(i,:);
        b(i+1) = b(i+1) - b(i);
    end
    
    % Back substitution
    R(m) = b(m)/A(m,m);
    for i = (m-1):-1:1
        sum = b(i);
        for j = (i+1):m
            sum = sum - A(i,j)*R(j);
        end
        R(i) = sum/A(i,i);
    end
end
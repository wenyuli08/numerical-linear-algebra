% Gauss Elimination with Scaled Partial Pivoting
function [R_GEwSPP] = GEwSPP(A,b)
    [m, n] = size(A);
    [rb, cb] = size(b);
    l = zeros(m,1);
    s = zeros(m,1);
    R_GEwSPP = zeros(m,1);
    
    % data preparation
    if m ~= n
        disp("A must be a square matrix.");
    end
    if rb ~= m || cb ~= 1
        disp("b must be a m x 1 matrix.");
    end
    d = det(A);
    if d == 0
        disp("Warning: singular matrices have either none or infinite number of solutions.")
    end

    % initialize the index array, li = i
    for i = 1:m
        l(i) = i;
        smax = 0;
        % calculate the scale array, si
        for j = 1:m
            smax = max(smax, abs(A(i,j)));
        end
        s(i) = smax;
    end
    
    % Gauss procedure
    % k is the index of the column in which new 0's are to be created
    for k = 1:(n-1)
        rmax = 0;
        % find the pivot for k-th column
        for i = k:m
            scaled_ratio = abs(A(l(i),k)/s(l(i)));
            if scaled_ratio > rmax
                rmax = scaled_ratio;
                pivot = i;
            end
        end
        l([k pivot]) = l([pivot k]);
        % calculate intermediate matrix when k-th column coefficient is 0
        for i = (k+1):m
            xmult = A(l(i),k)/A(l(k),k);
            A(l(i),k) = xmult;  % store the multiplier, note that its none-zero coefficient here doesn't matter
            % calculate (k+1)-th row after multiplication
            for j = (k+1):n
                A(l(i),j) = A(l(i),j) - xmult*A(l(k),j);
            end
            % bcur = bcur - xmult*bpvt
            b(l(i)) = b(l(i)) - A(l(i), k)*b(l(k));
        end
    end
    
    disp(A)
    disp(b)
    
    % Solve procedure
    R_GEwSPP(m) = b(l(m))/A(l(m),n);   % solve the equation w/ 1 unknown x
    for i = (n-1):-1:1
        sum = b(l(i));
        for j = (i+1):n
            sum = sum - A(l(i),j)*R_GEwSPP(j);
        end
        R_GEwSPP(i) = sum/A(l(i),i);
    end
end

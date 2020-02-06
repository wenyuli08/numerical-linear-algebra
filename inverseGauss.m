function [R] = inverseGauss(A)
    [m, n] = size(A);
    
    % data preparation
    if m ~= n
        disp("A must be a square matrix.");
    end
    d = det(A);
    if d == 0
        disp("Warning: singular matrices")
    end
    
    l = zeros(m,1); % row index vector
    s = zeros(m,1); % scale array
    R = zeros(n,n); % the result matrix
    I = eye(n); % the identity matrix

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
            I(l(i),:) = I(l(i),:) - A(l(i), k)*I(l(k),:);
        end
    end
   
    % Solve procedure
    for k = 1: n
        for i = n:-1:1
            for j = i:n
                I(l(i),k) = I(l(i),k) - A(l(i),j)*R(j,k);
            end
            R(i,k) = I(l(i),k)/A(l(i),i);
        end
    end
    
    disp("The inverse of A is: (using modified Gauss and Solve)")
end
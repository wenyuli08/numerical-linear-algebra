function [L,U] = tridiagonalB(A)
    [m, n] = size(A);
    if m ~= n || n ~= 4
        disp("A must be a 4x4 square matrix.")
    end
    
    % create empty arrays l(nxn) and u(nxn)
    l = zeros(n,n); % unit lower tridiagonal matrix, L
    u = zeros(n,n); % upper bidiagonal matrix, U
    
    % n = 4
    u(1,1) = A(1,1);
    for i = 2:4
        u(i-1,i) = A(i-1,i);
        l(i,i-1) = A(i,i-1)/u(i-1,i-1);
        u(i,i) = A(i,i) - l(i,i-1)*u(i-1,i);
    end
    % A = LU
    L = l;
    U = u;

    disp("L: ")
    disp(L)
    disp("U: ")
    disp(U)
end
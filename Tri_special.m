% Solve the tridiagonal system in which ai = ci = 1 for all i
function [R] = Tri_special(d,b)
    if length(d)~= length(b)
        disp("d and b must be of the same length.")
    end
    
    % a special case of procedure Tri: di = di - (1/di-1)*1; bi = bi -
    % (1/di-1)*bi-1
    n = length(d);
    for i = 2:n
        xmult = 1/d(i-1);
        d(i) = d(i) - xmult;
        b(i) = b(i) - xmult*b(i-1);
    end
    
    % Back substitution
    R(n) = b(n)/d(n);
    for i = (n-1):-1:1
        R(i) = (b(i) - R(i+1))/d(i);
    end
end
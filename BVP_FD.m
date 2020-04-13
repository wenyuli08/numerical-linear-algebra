% Finite Difference method for BVP
function BVP_FD(int,n)

    ta = int(1);
    tb = int(2);
    alpha = 1;
    beta = 0.5;

    u = @(t) 0;
    v = @(t) (1-t)/(1+t)^2;
    w = @(t) 1/(1+t)^2;

    h = (tb-ta)/n;
    format long;

    for i = 1:n-1
        t = ta+i*h;
        a(i) = -(1+(h/2)*w(t));
        d(i) = 2 + (h^2)*v(t);
        c(i) = -(1-(h/2)*w(t));
        b(i) = -h^2*u(t);
    end

    b(1) = b(1)-a(1)*alpha;
    b(n-1) = b(n-1) - c(n-1)*beta;

    a(n) = 0;
    for i=1:n-1
        a(i) = a(i+1);
    end
    [x] = Tri(a,d,c,b);
    [y] = [alpha x beta];
    tt = (ta:h:tb);

    plot(tt,y)

    error = exp(ta)-3*cos(ta)-alpha;
    disp('t-Value=')
    disp(ta)
    disp('Solution=')
    disp(exp(ta)-3*cos(ta))
    disp('Error=')
    disp(error)

    for i=9:9:n-1
        t = ta+i*h;
        error = exp(t)-3*cos(t)-x(i);
        disp('t-Value=')
        disp(t)
        disp('Solution=')
        disp(x(i))
        disp('Error=')
        disp(error)
    end

    error = exp(tb)-3*cos(tb)-beta;
    disp('t-Value=')
    disp(tb)
    disp('Solution=')
    disp(beta)
    disp('Error=')
    disp(error)
end
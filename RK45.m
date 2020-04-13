function [t,x,e] = RK45(h, f, t, x)
    c20 = 0.25;
    c21 = 0.25;
    c30 = 0.375;
    c31 = 0.09375;
    c32 = 0.28125;
    c40 = 12/13;
    c41 = 1932/2197;
    c42 = -7200/2197;
    c43 = 7296/2197;
    c51 = 439/216;
    c52 = -8;
    c53 = 3680/513;
    c54 = -845/4104;
    c60 = 0.5;
    c61 = -8/27;
    c62 = 2;
    c63 = -3544/2565;
    c64 = 1859/4104;
    c65 = -0.275;
    a1 = 25/216;
    a3 = 1408/2565;
    a4 = 2197/4104;
    a5 = -0.2;
    b1 = 16/135;
    b3 = 6656/12825;
    b4 = 28561/56430;
    b5 = -0.18;
    b6 = 2/55;
    K1 = h*f(t,x);
    K2 = h*f(t + c20*h, x + c21*K1);
    K3 = h*f(t + c30*h, x + c31*K1 + c32*K2);
    K4 = h*f(t + c40*h, x + c41*K1 + c42*K2 + c43*K3);
    K5 = h*f(t + h, x + c51*K1 + c52*K2 + c53*K3 + c54*K4);
    K6 = h*f(t + c60*h, x + c61*K1 + c62*K2 + c63*K3 + c64*K4 + c65*K5);
    x4 = x + a1*K1 + a3*K3 + a4*K4 + a5*K5;
    x = x + b1*K1 + b3*K3 + b4*K4 + b5*K5 + b6*K6;
    t = t + h;
    e = abs(x-x4);
end
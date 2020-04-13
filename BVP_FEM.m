function [ c ] = BVP_FEM(int,ya,yb,n)
%Finite Element solution if linear BVP
%Input: interval, boundary value, number of steps
%Output: solution values c
%Example usage: c = BVP_FEM([0 1],1,3,9)
a = int(1);
b = int(2);
h = (b-a)/(n+1);
alpha = (8/3)*h+2/h;
beta = (2/3)*h - 1/h;
format long
d = [-ya*beta;zeros(n-2,1);-yb*beta];
[cc] = Tri(beta*ones(n-1),alpha*ones(n),beta*ones(n-1),d);
c = [ya cc yb];
t = a:h:b;
disp('At the grid points')
disp(t)
yy = ((3-exp(-2))/(exp(2)-exp(-2)))*exp(2*t) + (exp(2)-3)/(exp(2)-exp(-2))*exp(-2*t); 
disp('the exact solution is')
disp(yy)
disp('and the FEM solution is')
disp(c)
hold all
plot(t,yy,'.','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);

plot(t,c+.03,'*','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
error = yy - c;
disp('with the error')
disp(error)
end

